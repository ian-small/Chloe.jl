

include("annotate_genomes.jl")
using JuliaWebAPI
using ArgParse
using Logging

const levels = Dict("info"=>Logging.Info, "debug"=> Logging.Debug, "warn" => Logging.Warn, 
"error"=>Logging.Error)

const ADDRESS = "tcp://127.0.0.1:9999"


MayBeIO = Union{Nothing,IOStream}


function chloe_svr(;refsdir = "reference_1116", address=ADDRESS,
    template = "optimised_templates.v2.tsv", level="warn", async=false, logfile::MayBeString=nothing)
    io::MayBeIO = nothing

    if logfile === nothing
        logger = ConsoleLogger(stderr,get(levels, level, Logging.Warn))
    else
        io = open(logfile::String, "a")
        logger = SimpleLogger(io, get(levels, level, Logging.Warn))
    end

    with_logger(logger) do
        reference = readReferences(refsdir, template)
        @info show_reference(reference)
        @info "using $(Threads.nthreads()) threads"

        function chloe(fasta::String, fname::MayBeString)
            annotate_one(fasta, reference, fname)
            return fname
        end

        function ping()
            return "OK"
        end
    
        process(
            JuliaWebAPI.create_responder([
                (chloe, false),
                (ping, false)

            ], address, true, "chloe"); async=async
        )

        if io !== nothing
            @info "closing logfile: $(logfile)"
            close(io::IOStream)
        end
    end
end

args = ArgParseSettings(prog="Chloë", autofix_names = true)  # turn "-" into "_" for arg names.

@add_arg_table! args begin
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
        default = ADDRESS
        help = "ZMQ address to listen on"
    "--logfile"
        arg_type=String
        metavar="FILE"
        help="log to file"
    "--level", "-l"
        arg_type = String
        metavar = "LOGLEVEL"
        default ="warn"
        help = "log level (warn,debug,info,error)"
        "--async"
        action = :store_true
        help = "run APIresponder async"
end
args.epilog = """
Run Chloe as a background ZMQ service
"""

function real_main() 
    parsed_args = parse_args(ARGS, args; as_symbols = true)
    # filter!(kv->kv.second ∉ (nothing, false), parsed_args)
    chloe_svr(;parsed_args...)
end


if abspath(PROGRAM_FILE) == @__FILE__
    real_main()
end

