


include("annotate_genomes.jl")
using JuliaWebAPI
using ArgParse
using Logging
using LogRoller
using Distributed

const LEVELS = Dict("info" => Logging.Info, "debug" => Logging.Debug, 
                    "warn" => Logging.Warn, "error"=> Logging.Error)

# const ADDRESS = "tcp://127.0.0.1:9999"
const ADDRESS = "ipc:///tmp/chloe-worker"

# change this if you change the API!
const VERSION = "1.0"

function git_version()
    try
        strip(read(`git rev-parse HEAD`, String))
    catch
        "unknown"
    end
end


MayBeString = Union{Nothing, String}
function chloe_distributed(;refsdir = "reference_1116", address=ADDRESS,
    template = "optimised_templates.v2.tsv", level="warn", nprocs=3,
    logfile::MayBeString=nothing)
    
    llevel = get(LEVELS, level, Logging.Warn)

    if logfile === nothing
        logger = ConsoleLogger(stderr,llevel)
    else
        logger = RollingLogger(logfile::String, 10 * 1000000, 2, llevel);
    end

    with_logger(logger) do

        machine = gethostname()
        reference = readReferences(refsdir, template)
        
        nthreads = Threads.nthreads()

        @info show_reference(reference)
        @info "chloe version $(VERSION) (git: $(git_version()[1:7])) threads=$(nthreads) on machine $(machine)"
        @info "connecting to $(address)"

        function chloe(fasta::String, fname::MayBeString)
            start = now()
            filename, target_id = fetch(@spawnat :any annotate_one(fasta, reference, fname))
            elapsed = now() - start
            @info "finished $(target_id) after $(elapsed)"
            return Dict("elapsed" => Dates.toms(elapsed), "filename" => filename, "ncid" => string(target_id))
        end

        function annotate(fasta::String)
            start = now()
            input = IOBuffer(fasta)
            io, target_id = fetch(@spawnat :any annotate_one(input, reference))
            sff = String(take!(io))
            elapsed = now() - start
            @info "finished $(target_id) after $(elapsed)"

            return Dict("elapsed" => Dates.toms(elapsed), "sff" => sff, "ncid" => string(target_id))
        end
        function ping()
            return "OK version=$(VERSION) procs=$(nprocs) on $(machine)"
        end
        function threads()
            return nprocs
        end
        # we need to create separate ZMQ sockets to ensure strict
        # request/response (not e.g. request-request response-response)
        @sync for i in 1:nprocs
            @async process(
                JuliaWebAPI.create_responder([
                        (chloe, false),
                        (annotate, false),
                        (ping, false),
                        (threads, false)

                    ], address, false, "chloe"); async=false
                )

        end

        if isa(logger, RollingLogger)
            close(logger::RollingLogger)
        end

    end
end

distributed_args = ArgParseSettings(prog="Chloë", autofix_names = true)  # turn "-" into "_" for arg names.

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
        metavar ="URL"
        help = "ZMQ DEALER address to connect to"
    "--logfile"
        arg_type=String
        metavar="FILE"
        help="log to file"
    "--level", "-l"
        arg_type = String
        metavar = "LOGLEVEL"
        default ="warn"
        help = "log level (warn,debug,info,error)"
    "--nprocs"
        arg_type = Int
        default = 3
        help = "number of distributed processes"

end
distributed_args.epilog = """
Run Chloe as a background ZMQ service with distributed annotation processes
"""

distributed_args = parse_args(ARGS, distributed_args; as_symbols = true)

procs = addprocs(get(distributed_args, :nprocs, 3); topology=:master_worker)
# delete!(distributed_args,:nprocs)

@info "processes: $(procs)"
Sys.set_process_title("chloe-distributed")

# sic! src/....
@everywhere procs include("src/annotate_genomes.jl")

chloe_distributed(;distributed_args...)




