module ChloeDistributed
export distributed_main, annotate_one_task

import Distributed
import Distributed: addprocs, rmprocs, @spawnat, @everywhere, nworkers


import Base
import JuliaWebAPI: APIResponder, APIInvoker, apicall, ZMQTransport, register, process
import ArgParse: ArgParseSettings, @add_arg_table!, parse_args
import Dates: now, toms
import Crayons: @crayon_str
import StringEncodings: encode
import ZMQ

import ..WebAPI: TerminatingJSONMsgFormat
import ..Annotator: MayBeString, verify_refs, ChloeConfig, default_gsrefsdir
import ..ZMQLogging: set_global_logger
import ..Broker: check_endpoints, remove_endpoints

include("../globals.jl")
include("dist_globals.jl")
include("tasks.jl")
function git_version()
    repo_dir = dirname(@__FILE__)
    try
        # older version of git don't have -C
        # strip(read(pipeline(`git -C "$repo_dir" rev-parse HEAD`, stderr=devnull), String))
        strip(read(pipeline(`sh -c 'cd "$repo_dir" && git rev-parse HEAD'`, stderr=devnull), String))
    catch e
        "unknown"
    end
end
# from https://stackoverflow.com/questions/27677399/julia-how-to-copy-data-to-another-processor-in-julia
# https://github.com/ChrisRackauckas/ParallelDataTransfer.jl
# function sendto(p::Int; args...)
#     for (nm, val) in args
#         @spawnat(p, Base.eval(Main, Expr(:(=), nm, val)))
#     end
# end

# function sendto(p::Int; args...)

#     function send(nm, val)
#         @spawnat p Base.eval(Main, Expr(:(=), nm, val))
#     end

#     [send(nm, val) for (nm, val) in args]

# end

# getfrom(p::Int, nm::Symbol; mod=Main) = fetch(@spawnat(p, getfield(mod, nm)))


# function sendto(ps::Vector{Int}; args...)
#     for p in ps
#         sendto(p; args...)
#     end
# end

# function sendto(ps::Vector{Int}; args...)
#     [tsk for p in ps for tsk in sendto(p; args...) ]
# end

function create_responder(apispecs::Vector{Function}, addr::String, ctx::ZMQ.Context)
    api = APIResponder(ZMQTransport(addr, ZMQ.REP, false, ctx), TerminatingJSONMsgFormat(), "chloe", false)
    for func in apispecs
        register(api, func)
    end
    api
end

function arm_procs_full(procs, backend::MayBeString=nothing, level::String="info";
    gsrefsdir="default", template="default")

    @everywhere procs begin
        # If we have a raw 'import Chloe' here then
        # precompilation (of the Chloe package) tries to recurse and compile Chloe *again* (I think) and fails.
        # "hiding" the import inside a quote seems to work.
        eval(quote
            import Chloe
        end)
        Chloe.set_global_logger($level, $backend; topic="annotator")
        global REFERENCE = Chloe.ReferenceDb(; gsrefsdir=$gsrefsdir, template=$template)
    end

end

function chloe_distributed(; gsrefsdir="default", address=ZMQ_WORKER,
    template="default", level="warn", workers=3,
    backend::MayBeString=nothing, broker::MayBeString=nothing, reference_dir::MayBeString=nothing)

    if !isnothing(backend)
        if backend == "default"
            backend = ZMQ_BACKEND
        end
    end
    set_global_logger(level, backend; topic="annotator")
    if !isnothing(reference_dir)
        gsrefsdir = normpath(joinpath(expanduser(reference_dir), "gsrefs"))
        template = normpath(joinpath(dirname(gsrefsdir), DEFAULT_TEMPLATE))
    else
        if gsrefsdir == "default"
            gsrefsdir = default_gsrefsdir()
        end
        if template == "default"
            template = normpath(joinpath(dirname(gsrefsdir), DEFAULT_TEMPLATE))
        end
    end

    # don't wait for workers to find the wrong directory
    verify_refs(gsrefsdir, template)

    # user may have added run with 
    # julia command -p2 etc.
    # don't send to process 1
    procs = filter(w -> w != 1, Distributed.workers())
    toadd = workers - length(procs)
    if toadd > 0
        addprocs(toadd; topology=:master_worker, exeflags="--project=$(Base.active_project())")
    end

    procs = filter(w -> w != 1, Distributed.workers())

    function arm(new_procs)
        arm_procs_full(new_procs, backend, level; gsrefsdir=gsrefsdir, template=template)
    end

    arm(procs)

    chloe_listen(address, broker, arm)
end

function chloe_listen(address::String, broker::MayBeString, arm_procs::Function)
    success = crayon"bold green"
    procs = Distributed.workers()

    nlisteners = length(procs)


    pid = getpid()
    machine = gethostname()
    git = git_version()[1:7]
    nthreads = Threads.nthreads()

    @info "starting: chloe version $VERSION (git: $git) pid=$pid workers=$nlisteners threads=$nthreads on machine $machine"
    @info "connecting to $address"
    # @info reference

    # All sent to nlisteners I don't need this anymore (except see add_nlisteners below)
    # reference = nothing # force garbage collection
    # GC.gc()
    nannotations = 0

    function chloe(fasta::String, outputsff::MayBeString=nothing, config::Union{Nothing,Dict{String,V} where V<:Any}=nothing, task_id::MayBeString=nothing)
        start = now()
        cfg = if isnothing(config)
            ChloeConfig()
        else
            ChloeConfig(config)
        end
        filename, target_id = fetch(@spawnat :any Main.Chloe.annotate_one_task(fasta, outputsff, task_id, cfg))
        elapsed = now() - start
        @info success("finished $target_id after $elapsed")
        nannotations += 1
        return Dict("elapsed" => toms(elapsed), "filename" => filename, "ncid" => string(target_id), "config" => cfg)
    end

    function batch_annotate(directory::String, config::Union{Nothing,Dict{String,V} where V<:Any}=nothing, task_id::MayBeString=nothing)
        start = now()
        cfg = if isnothing(config)
            ChloeConfig()
        else
            ChloeConfig(config)
        end
        n = fetch(@spawnat :any annotate_batch_task(directory, task_id, cfg))
        elapsed = now() - start
        @info success("finished $directory after $elapsed")
        nannotations += n
        return Dict("elapsed" => toms(elapsed), "directory" => directory, "config" => cfg)
    end

    function decompress(fasta::String)
        # decode a latin1 encoded binary string
        read(encode(fasta, "latin1") |> IOBuffer |> GzipDecompressorStream, String)
    end

    function annotate(fasta::String, config::Union{Nothing,Dict{String,V} where V<:Any}=nothing, task_id::MayBeString=nothing)
        start = now()
        if startswith(fasta, "\u1f\u8b")
            # assume latin1 encoded binary gzip file
            n = length(fasta)
            fasta = decompress(fasta)
            @debug "decompressed fasta length $(n) -> $(length(fasta))"
        end
        cfg = if isnothing(config)
            ChloeConfig()
        else
            ChloeConfig(config)
        end
        input = IOBuffer(fasta)
        io, target_id = fetch(@spawnat :any Main.Chloe.annotate_one_task(input, task_id, cfg))
        sff = String(take!(io))
        elapsed = now() - start
        @info success("finished $target_id after $elapsed")
        nannotations += 1

        return Dict("elapsed" => toms(elapsed), "sff" => sff, "ncid" => string(target_id), "config" => cfg)
    end

    function ping()
        return "OK version=$VERSION git=$git #anno=$nannotations pid=$pid threads=$nthreads workers=$nlisteners on $machine"
    end

    # `bin/chloe.py terminate` uses this to find out how many calls of :terminate
    # need to be made to stop all responders. It's hard to cleanly
    # stop process(APIResponder) from the outside since it is block wait on 
    # the zmq sockets.
    function nconn()
        return nlisteners
    end

    function bgexit(endpoint::String, n::Int)
        i = APIInvoker(endpoint)
        maxfails = 10
        while n > 0 && maxfails > 0
            # main task has removed nlisteners
            if nlisteners === 0
                break
            end

            # the real problem here is:
            # -- unless we are running our own broker --
            # there is no way to ensure we are terminating
            # *our* listeners... we send a pid so it can check

            # @async to allow for main task to count down nlisteners
            res = fetch(@async apicall(i, ":exit", pid))
            code = res["code"]
            @debug "code=$(code) workers=$(nlisteners)"
            # ping wrong server... oops
            if code === 200
                n -= 1
            else
                maxfails -= 1
            end
        end
    end
    function exit(endpoint::MayBeString=nothing)
        # use broker url if any
        if endpoint === nothing
            endpoint = broker
        end
        if endpoint === nothing
            error("No broker endpoint!")
        end
        @async bgexit(endpoint, nlisteners)
        return "Exit scheduled for $(nlisteners) workers"
    end

    function add_workers(n::Int, endpoint::MayBeString=nothing)
        if n < 0
            if -n >= length(procs)
                error("use 'exit' to exit chloe!")
            end
            # remove listeners & workers
            @async begin
                m = length(procs) + n
                remove = procs[m+1:end]
                rmprocs(remove) # wait
                # update globals !not nlisteners thou.
                procs = procs[1:m]
                # send terminate to abs(n) listeners
                if endpoint === nothing
                    endpoint = broker
                end
                if endpoint !== nothing
                    # main task will update nlisteners count
                    # when they exit
                    bgexit(endpoint, abs(n))
                end
                @info "removed $(remove) processes using endpoint $(endpoint)"
            end

        elseif n > 0
            # add workers
            @async begin
                # ensure topology is the same
                added = addprocs(n, topology=:master_worker, exeflags="--project=$(Base.active_project())")
                arm_procs(added)
                @info "added $(added) processes"
                # update globals
                procs = [procs..., added...]
                nlisteners = length(procs)
                for w in added
                    @async bgworker(w)
                end
            end
        end
        msg = n > 0 ? "add" : "remove"
        "Chloë is scheduled to $(msg) $(abs(n)) worker$(abs(n) !== 1 ? "s" : "")"
    end
    # we need to create separate ZMQ sockets to ensure strict
    # request/response (not e.g. request-request response-response)
    # we expect to *connect* to a ZMQ DEALER/ROUTER (see bin/broker.py
    # or src/broker.jl) that forms the actual front end.

    ctx = ZMQ.Context()

    function cleanup()
        # zmq listeners are already finished
        # only worker processes need destroying
        close(ctx)
        try
            rmprocs(procs, waitfor=20)
        catch e
            @warn "background processes took too long to exit $e"
        end
    end

    atexit(cleanup)

    done = Channel{Int}()

    function bgworker(workno::Int)
        @debug "starting worker $workno"
        process(
            create_responder([
                    chloe,
                    annotate,
                    # batch_annotate,
                    ping,
                    nconn,
                    exit,
                    add_workers,
                ], address, ctx)
        )
        # :terminate called so process loop is finished
        put!(done, workno)
    end

    # kick off worker tasks to listen
    # on zmq endpoints
    for workno in procs
        @async bgworker(workno)
    end
    @info "Chloë is ready to accept requests 🍰!"
    while nlisteners > 0
        w = take!(done) # block here until listeners exit (by putting proc num on queue)
        nlisteners -= 1
        @debug "worker done: $w workers=$nlisteners"
    end

    @info success("Chloë done: pid=$(pid) exiting.....🍹")
    # atexit cleanup and kill(broker) task should be called now

end

function get_distributed_args(args::Vector{String}=ARGS)
    distributed_args = ArgParseSettings(prog="Chloë", autofix_names=true)  # turn "-" into "_" for arg names.

    @add_arg_table! distributed_args begin
        "--reference", "-r"
        arg_type = String
        dest_name = "reference_dir"
        metavar = "DIRECTORY"
        help = "reference directory (takes precedence over --gsrefs and --template options)"
        "--gsrefs", "-g"
        arg_type = String
        default = "default"
        dest_name = "gsrefsdir"
        metavar = "DIRECTORY"
        help = "gsrefs reference directory [default: $(DEFAULT_GSREFS)]"
        "--template", "-t"
        arg_type = String
        default = "default"
        metavar = "TSV"
        dest_name = "template"
        help = "template tsv [default: {reference directory}/$(DEFAULT_TEMPLATE)]"
        "--address", "-a"
        arg_type = String
        metavar = "URL"
        default = ZMQ_WORKER
        help = "ZMQ DEALER address to connect to"
        "--level", "-l"
        arg_type = String
        metavar = "LOGLEVEL"
        default = "info"
        help = "log level (warn,debug,info,error,critical)"
        "--workers", "-w"
        arg_type = Int
        default = 3
        help = "number of distributed processes"
        "--broker", "-b"
        arg_type = String
        metavar = "URL"
        help = "run a broker in the background (use \"default\" for default endpoint)"
        "--backend", "-z"
        arg_type = String
        metavar = "URL"
        help = "log to zmq endpoint"


    end

    distributed_args.epilog = """
    Run Chloë as a ZMQ service with distributed annotation processes.
    Requires a ZMQ DEALER/ROUTER to connect to unless `--broker, -b` specifies
    an endpoint in which case it runs its own broker.

    You can also use julia arguments to specify workers with e.g.
    `julia -p4 ....` etc.
    """
    parse_args(args, distributed_args; as_symbols=true)

end

function run_broker(worker::String=ZMQ_WORKER, client::String=ZMQ_ENDPOINT)
    #  see https://discourse.julialang.org/t/how-to-run-a-process-in-background-but-still-close-it-on-exit/27231
    # src = dirname(@__FILE__)
    julia = joinpath(Sys.BINDIR, "julia")
    if !Sys.isexecutable(julia)
        error("Can't find julia executable to run broker, best guess: $julia")
    end
    msg = check_endpoints(worker, client)
    if msg !== nothing
        @error "$(msg). Only need one broker running"
        Base.exit(1) # can we throw....
    end
    # cmd = `$julia --project=$(pwd()) -q --startup-file=no "$src/broker.jl" --worker=$worker --client=$client`
    cmd = `$julia --project="$(Base.active_project())" -q --startup-file=no -e "import Chloe; Chloe.broker_main()"  -- --worker=$worker --client=$client`
    # @info "running broker as: $cmd"
    # wait = false means stdout,stderr are connected to /dev/null
    task = run(cmd; wait=false)
    atexit(() -> kill(task))
    remove_endpoints(worker, client)
    task
end


function maybe_launch_broker(distributed_args)
    function find_endpoint()
        endpoint = tmplt = "/tmp/chloe-worker"
        n = 0
        while isfile(endpoint)
            n += 1
            endpoint = "$(tmplt)$(n)"
        end
        "ipc://$(endpoint)"
    end

    client_url = get(distributed_args, :broker, nothing)

    if client_url !== nothing
        if client_url == "default"
            client_url = ZMQ_ENDPOINT
            distributed_args[:broker] = client_url
        end
        if startswith(client_url, "@")
            # hack just so the server knows how to terminate itself
            distributed_args[:broker] = client_url[2:end]
        else
            # really start a broker
            if distributed_args[:address] === client_url
                distributed_args[:address] = find_endpoint()
                @warn "broker and worker endpoints clash: redirecting workers to $(distributed_args[:address])"
            end

            @info "Starting broker. Connect to: $client_url"
            run_broker(distributed_args[:address], client_url)
        end
    end
    distributed_args
end

function distributed_main(args::Vector{String}=ARGS)
    Sys.set_process_title("chloe-distributed")
    distributed_args = get_distributed_args(args)
    distributed_args = maybe_launch_broker(distributed_args)

    chloe_distributed(; distributed_args...)
end

end # module
