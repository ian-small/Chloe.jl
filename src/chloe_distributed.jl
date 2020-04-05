include("annotate_genomes.jl")
include("ZMQLogger.jl")
using JuliaWebAPI
using ArgParse
# using LogRoller
using Distributed
using Crayons
import StringEncodings: encode

const success = crayon"bold green"
const ADDRESS = "tcp://127.0.0.1:9999"
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
        if !startswith(fasta, '>')
            # assume latin1 encoded binary
            @info "compressed fasta length $(length(fasta))"
            fasta = read(GzipDecompressorStream(IOBuffer(encode(fasta, "latin1"))), String)
            @info "decompressed fasta length $(length(fasta))"
        end
        start = now()

        input = IOContext(IOBuffer(fasta), stdin)

        io, target_id = fetch(@spawnat :any annotate_one(reference, input))
        sff = String(take!(io))
        elapsed = now() - start
        @info success("finished $target_id after $elapsed")

        return Dict("elapsed" => Dates.toms(elapsed), "sff" => sff, "ncid" => string(target_id))
    end

    function ping()
        return "OK version=$VERSION git=$git threads=$nthreads procs=$nprocs on $machine"
    end

    function nconn()
        return nprocs
    end
    # we need to create separate ZMQ sockets to ensure strict
    # request/response (not e.g. request-request response-response)
    # we expect to *connect* to a ZMQ DEALER/ROUTER (see bin/broker.py)
    # that forms the actual front end.
    @sync for i in 1:nprocs
        @async process(
            JuliaWebAPI.create_responder([
                    (chloe, false),
                    (annotate, false),
                    (ping, false),
                    (nconn, false)

                ], address, false, "chloe"); async = false
            )

    end

end

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
    "--logendpoint"
        arg_type = String
        metavar = "ZMQ"
        help = "log to zmq endpoint"
    "--level", "-l"
        arg_type = String
        metavar = "LOGLEVEL"
        default = "info"
        help = "log level (warn,debug,info,error)"
    "--nprocs"
        arg_type = Int
        default = 3
        help = "number of distributed processes"

end
distributed_args.epilog = """
Run Chloe as a background ZMQ service with distributed annotation processes.
Requires a ZMQ DEALER/ROUTER to connect to.
"""

distributed_args = parse_args(ARGS, distributed_args; as_symbols = true)

# delete!(distributed_args,:nprocs)

Sys.set_process_title("chloe-distributed")

chloe_distributed(;distributed_args...)
