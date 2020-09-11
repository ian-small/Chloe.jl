

include("annotate_genomes.jl")
include("sff2GFF3.jl")
include("make_suffix_array.jl")
import ArgParse: ArgParseSettings, @add_arg_table!, parse_args
import Logging

const levels = Dict("info" => Logging.Info, "debug" => Logging.Debug, "warn" => Logging.Warn, 
"error" => Logging.Error)
 
function chloe(;refsdir="reference_1116", fasta_files=String[],
    template="optimised_templates.v2.tsv", output::MayBeString=nothing)

    annotate(refsdir, template, fasta_files, output)

end


# const julia_v07 = VERSION > v"0.7-"

cmd_args = ArgParseSettings(prog="Chloë", autofix_names=true)  # turn "-" into "_" for arg names.

@add_arg_table! cmd_args begin
    "annotate"
        help = "annotate fasta files"
        action = :command
    "gff3"
        help = "write out gff3 files from sff files"
        action = :command
    "suffix"
        help = "generate suffix arrays from fasta files"
        action = :command
    # "jld"
    #     help = "generate jld reference file"
    #     action = :command
    "--level", "-l"
        arg_type = String
        default = "warn"
        help = "log level (warn,debug,info,error)"
end

@add_arg_table! cmd_args["gff3"]  begin
    "sff-files"
        arg_type = String
        nargs = '+'
        required = true
        action = :store_arg
        help = "sff files to process" 
        "--directory", "-d"
        arg_type = String
        help = "output directory" 
end

@add_arg_table! cmd_args["suffix"]  begin
    "fasta-files"
        arg_type = String
        nargs = '+'
        required = true
        action = :store_arg
        help = "fasta files to process"
end

@add_arg_table! cmd_args["annotate"]  begin
    "fasta-files"
        arg_type = String
        nargs = '+'
        required = true
        action = :store_arg
        help = "fasta files to process"
    "--output", "-o"
        arg_type = String
        help = "output filename (or directory if multiple fasta files)"
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
end

# @add_arg_table! cmd_args["jld"]  begin
#     "--output", "-o"
#         arg_type = String
#         help = "output filename"
#         required = true
#     "--reference", "-r"
#         arg_type = String
#         default = "reference_1116"
#         dest_name = "refsdir"
#         metavar = "DIRECTORY"
#         help = "reference directory"
#     "--template", "-t"
#         arg_type = String
#         default = "optimised_templates.v2.tsv"
#         metavar = "TSV"
#         dest_name = "template"
#         help = "template tsv"
# end

# args.epilog = """
#     examples:\n
#     \ua0\ua0 # chloe.jl -t template.tsv -r reference_dir fasta1 fasta2 ...\n
#     """

function cmd_main() 
    parsed_args = parse_args(ARGS, cmd_args; as_symbols=true)
    level = parsed_args[:level]
    Logging.with_logger(Logging.ConsoleLogger(stderr, get(levels, level, Logging.Warn))) do 
        cmd = parsed_args[:_COMMAND_]
        a = parsed_args[cmd]
        if cmd == :gff3
            writeallGFF3(;a...)
        elseif cmd == :annotate
            chloe(;a...)
        elseif cmd == :suffix
            writesuffixarray(;a...)
        end
    end

end


if abspath(PROGRAM_FILE) == @__FILE__
    cmd_main()
end

