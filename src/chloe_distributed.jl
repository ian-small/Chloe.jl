include("annotate_genomes.jl")
include("ZMQLogger.jl")
# include("broker.jl")
using JuliaWebAPI
using ArgParse
# using LogRoller
using Distributed
using Crayons
import StringEncodings: encode

const success = crayon"bold green"
const ADDRESS = "tcp://127.0.0.1:9467"
# const ADDRESS = "ipc:///tmp/chloe-worker"

# change this if you change the API!
const VERSION = "1.0"

function git_version()
    try
        strip(read(`git rev-parse HEAD`, String))
    catch
        "unknown"
    end
end

# from 
function exit_on_sigint(on::Bool)
    # from https://github.com/JuliaLang/julia/pull/29383
    # and https://github.com/JuliaLang/julia/pull/29411
    ccall(:jl_exit_on_sigint, Cvoid, (Cint,), on)
end

function create_responder(apispecs::Array{Function}, addr::String, ctx::Context)
    api = APIResponder(ZMQTransport(addr, REP, false, ctx), JSONMsgFormat(), "chloe", false)
    for func in apispecs
        register(api, func)
    end
    api
end

function chloe_distributed(;refsdir = "reference_1116", address = ADDRESS,
    template = "optimised_templates.v2.tsv", level = "warn", nprocs = 3,
    logendpoint::MayBeString = nothing)

    procs = addprocs(nprocs; topology = :master_worker)
    # sic! src/....
    @everywhere procs include("src/annotate_genomes.jl")
    @everywhere procs include("src/ZMQLogger.jl")
    # can't use rolling logger for procs because of file contentsion
    for p in procs
        @spawnat p set_global_logger(logendpoint, level; topic = "annotator")
    end
    set_global_logger(logendpoint, level; topic = "annotator")
    
    machine = gethostname()
    reference = readReferences(refsdir, template)
    git = git_version()[1:7]
    
    nthreads = Threads.nthreads()
    @info "processes: $nprocs"
    @info reference
    @info "chloe version $VERSION (git: $git) threads=$nthreads on machine $machine"
    @info "connecting to $address"

    function chloe(fasta::String, fname::MayBeString)
        start = now()
        filename, target_id = fetch(@spawnat :any annotate_one(reference, fasta, fname))
        elapsed = now() - start
        @info success("finished $target_id after $elapsed")
        return Dict("elapsed" => Dates.toms(elapsed), "filename" => filename, "ncid" => string(target_id))
    end

    function annotate(fasta::String)
        start = now()

        if !startswith(fasta, '>')
            # assume latin1 encoded binary
            @info "compressed fasta length $(length(fasta))"
            fasta = read(encode(fasta, "latin1") |> IOBuffer |> GzipDecompressorStream, String)
            # fasta = read(GzipDecompressorStream(IOBuffer(encode(fasta, "latin1"))), String)
            @info "decompressed fasta length $(length(fasta))"
        end

        input = IOContext(IOBuffer(fasta))

        io, target_id = fetch(@spawnat :any annotate_one(reference, input))
        sff = String(take!(io))
        elapsed = now() - start
        @info success("finished $target_id after $elapsed")

        return Dict("elapsed" => Dates.toms(elapsed), "sff" => sff, "ncid" => string(target_id))
    end

    function ping()
        return "OK version=$VERSION git=$git threads=$nthreads procs=$nprocs on $machine"
    end

    # `bin/chloe.py terminate` uses this to find out how many calls of :terminate
    # need to be made to stop all responders. It's hard to cleanly
    # stop process(APIResponder) from the outside since it is block wait on 
    # the zmq sockets.
    function nconn()
        return nprocs
    end

    # we need to create separate ZMQ sockets to ensure strict
    # request/response (not e.g. request-request response-response)
    # we expect to *connect* to a ZMQ DEALER/ROUTER (see bin/broker.py)
    # that forms the actual front end.
    ctx = Context()

    function cleanup()
        close(ctx)
        try
            rmprocs(procs, waitfor = 20)
        catch
        end
    end
    
    atexit(cleanup)

    @sync for p in procs
        @async JuliaWebAPI.process(
            create_responder([
                    chloe,
                    annotate,
                    ping,
                    nconn,
                ], address, ctx)
            )

    end
end

function args()
    distributed_args = ArgParseSettings(prog = "Chloë", autofix_names = true)  # turn "-" into "_" for arg names.

    @add_arg_table! distributed_args begin
        "--reference", "-r"
        arg_type = String
        default = "reference_1116"
        dest_name = "refsdir"
        metavar = "DIRECTORY"
        help = "reference directory"
        "--template", "-t"
        arg_type = String
        default = "optimised_templates.v2.tsv"
        metavar = "TSV"
        dest_name = "template"
        help = "template tsv"
        "--address", "-a"
        arg_type = String
        metavar = "URL"
        default = ADDRESS
        help = "ZMQ DEALER address to connect to"
        "--level", "-l"
        arg_type = String
        metavar = "LOGLEVEL"
        default = "info"
        help = "log level (warn,debug,info,error)"
        "--nprocs", "-n"
        arg_type = Int
        default = 3
        help = "number of distributed processes"
        "--broker"
            arg_type = String
            metavar = "URL"
            help = "run the broker in the background"
        "--logendpoint"
            arg_type = String
            metavar = "ZMQ"
            help = "log to zmq endpoint"
    end
    distributed_args.epilog = """
    Run Chloe as a ZMQ service with distributed annotation processes.
    Requires a ZMQ DEALER/ROUTER to connect to unless `--broker` specifies
    an endpoint in which case it runs its own broker.
    """
    parse_args(ARGS, distributed_args; as_symbols = true)

end

function run_broker(worker, client)
    #  see https://discourse.julialang.org/t/how-to-run-a-process-in-background-but-still-close-it-on-exit/27231
    src = dirname(@__FILE__)
    julia = joinpath(Sys.BINDIR, "julia")
    if !Sys.isexecutable(julia)
        error("Can't find julia executable to run broker, best guess: $julia")
    end
    cmd = `$julia -q --startup-file=no "$src/broker.jl" --worker=$worker --client=$client`
    # wait = false means stdout,stderr are connected to /dev/null
    task = run(cmd; wait = false)
    atexit(()->kill(task))
    task
    # open(pipeline(cmd))
end

function run_broker2(worker, client)
    # ugh! `@spawnat :any annotate...` will block on this process... which
    # will never return.
    procs = addprocs(1; topology = :master_worker)
    @everywhere procs include("src/broker.jl")
    @async fetch(@spawnat procs[1] run_broker(worker, client))
end
    
function main()
    # exit_on_sigint(false)
    Sys.set_process_title("chloe-distributed")
    distributed_args = args()
    client_url = pop!(distributed_args, :broker, nothing)


    if client_url !== nothing
        @info "Starting broker. Connect to: $client_url"
        run_broker(distributed_args[:address], client_url)
    end
    chloe_distributed(;distributed_args...)

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
