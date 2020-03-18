

include("annotate_genomes.jl")
using ArgParse
using Logging

const levels = Dict("info"=>Logging.Info, "debug"=> Logging.Debug, "warn" => Logging.Warn, 
"error"=>Logging.Error)
 
function chloe(;refsdir = "reference_1116", fasta_files = String[], verbose = false,
    template = "optimised_templates.v2.tsv", level="warn", output::MayBeString = nothing)

    with_logger(ConsoleLogger(stderr,get(levels, level, Logging.Warn))) do
        annotate(refsdir, template, fasta_files, output)
    end
end


# const julia_v07 = VERSION > v"0.7-"

args = ArgParseSettings(autofix_names = true)  # turn "-" into "_" for arg names.

@add_arg_table! args begin
    "fasta-files"
        arg_type = String
        nargs = '+'
        required = true
        action = :store_arg
        help = "fasta files to process"
    "--output", "-o"
        arg_type = String
        help = "output filename"
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
    "--verbose", "-v"
        action = :store_true
        help = "increase verbosity"
        "--level", "-l"
        arg_type = String
        default ="warn"
        help = "log level (warn,debug,info,error)"
end
args.epilog = """
    examples:\n
    \ua0\ua0 # chloe.jl -t template.tsv -r reference_dir fasta1 fasta2 ...\n
    """

function real_main() 
    parsed_args = parse_args(ARGS, args; as_symbols = true)
    # filter!(kv->kv.second ∉ (nothing, false), parsed_args)
    chloe(;parsed_args...)
end


if abspath(PROGRAM_FILE) == @__FILE__
    real_main()
end

