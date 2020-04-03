

include("annotate_genomes.jl")
using JuliaWebAPI
using ArgParse
using Logging
using LogRoller
using Dates

const LEVELS = Dict("info"=>Logging.Info, "debug"=> Logging.Debug, 
                    "warn" => Logging.Warn, "error"=>Logging.Error)

const ADDRESS = "tcp://127.0.0.1:9999"

# change this if you change the API!
const VERSION = "1.0"

function git_version()
    try
        strip(read(`git rev-parse HEAD`, String))
    catch
        "unknown"
    end
end

function chloe_svr(;refsdir = "reference_1116", address=[ADDRESS],
    template = "optimised_templates.v2.tsv", level="warn",
    logfile::MayBeString=nothing, connect=false, nconn=1)
    async = false
    
    llevel = get(LEVELS, level, Logging.Warn)
    
    if length(address) === 0
        push!(address, ADDRESS)
    end

    address = repeat(address, nconn)

    if logfile === nothing
        logger = ConsoleLogger(stderr,llevel)
    else
        logger = RollingLogger(logfile::String, 10 * 1000000, 2, llevel);
    end

    conn = connect ? "connecting to" : "listening on"

    with_logger(logger) do
        Sys.set_process_title("chloe-svr")
        machine = gethostname()
        reference = readReferences(refsdir, template)
        @info reference
        @info "chloe version $(VERSION) (git: $(git_version()[1:7])) using $(nconn) threads on machine $(machine)"
        @info "$(conn) $(address)"

        function chloe(fasta::String, fname::MayBeString)
            @info "running on thread: $(Threads.threadid())"
            start = now()
            filename, target_id = annotate_one(reference, fasta, fname)
            elapsed = now() - start
            @info "finished $(target_id) on thread: $(Threads.threadid()) after $(elapsed)"
            return Dict("elapsed" => Dates.toms(elapsed), "filename" => filename, "ncid" => string(target_id))
        end

        function annotate(fasta::String)
            @info "running on thread: $(Threads.threadid())"
            start = now()
            input = IOBuffer(fasta)
            io, target_id = annotate_one(reference, input)
            sff = String(take!(io))
            elapsed = now() - start
            @info "finished $(target_id) on thread: $(Threads.threadid()) after $(elapsed)"

            return Dict("elapsed" => Dates.toms(elapsed), "sff" => sff, "ncid" => string(target_id))
        end
        function ping()
            return "OK version=$(VERSION) on thread=$(Threads.threadid())/$(length(address)) on $(machine)"
        end
        function nconn()
            return length(address)
        end
        if length(address) == 1
            process(
                    JuliaWebAPI.create_responder([
                            (chloe, false),
                            (annotate, false),
                            (ping, false),
                            (nconn, false)

                        ], address[1], !connect, "chloe"); async=async
                    )
        else
            Threads.@threads for addr in address
                process(
                    JuliaWebAPI.create_responder([
                        (chloe, false),
                        (annotate, false),
                        (ping, false),
                        (nconn, false)

                    ], addr, !connect, "chloe"); async=async
                )
            end
        end
        if isa(logger, RollingLogger)
            close(logger::RollingLogger)
        end

    end
end

svr_args = ArgParseSettings(prog="Chloë", autofix_names = true)  # turn "-" into "_" for arg names.

@add_arg_table! svr_args begin
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
        default = []
        help = "ZMQ address(es) to listen on or connect to"
        action = :append_arg
    "--logfile"
        arg_type=String
        metavar="FILE"
        help="log to file"
    "--level", "-l"
        arg_type = String
        metavar = "LOGLEVEL"
        default ="warn"
        help = "log level (warn,debug,info,error)"
    "--connect", "-c"
        action = :store_true
        help = "connect to addresses instead of bind"
    "--nconn"
        arg_type = Int
        default = 1
        help = "number of threads when connecting"

end
svr_args.epilog = """
Run Chloe as a background ZMQ service
"""

function svr_main() 
    parsed_args = parse_args(ARGS, svr_args; as_symbols = true)
    # filter!(kv->kv.second ∉ (nothing, false), parsed_args)
    chloe_svr(;parsed_args...)
end


if abspath(PROGRAM_FILE) == @__FILE__
    svr_main()
end

