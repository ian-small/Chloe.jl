
module CmdLine
export cmd_main

import ArgParse: ArgParseSettings, @add_arg_table!, parse_args

include("globals.jl")

import ..Annotator
import ..Sff2Gff
import ..SuffixArray

 
function chloe(;refsdir=DEFAULT_REFS, fasta_files=String[],
    template=DEFAULT_TEMPLATE, output::Union{Nothing,String}=nothing)
    if refsdir == "default"
        refsdir = joinpath(HERE, "..", DEFAULT_REFS)
    end
    if template == "default"
        template = joinpath(HERE, "..", DEFAULT_TEMPLATE)
    end
    Annotator.annotate(refsdir, template, fasta_files, output)

end


function getargs()
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
        "mmap"
            help = "generate mmap files from gwsas/fasta files"
            action = :command
        "--level", "-l"
            arg_type = String
            default = "info"
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
        "--directory", "-d"
            arg_type = String
            help = "output directory" 
    end

    @add_arg_table! cmd_args["mmap"]  begin
        "gwsas_fasta"
            arg_type = String
            nargs = '+'
            required = true
            action = :store_arg
            help = "gwsas/fasta files to process"

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
            default = "default"
            dest_name = "refsdir"
            metavar = "DIRECTORY"
            help = "reference directory [default: $(DEFAULT_REFS)]"
        "--template", "-t"
            arg_type = String
            default = "default"
            metavar = "TSV"
            dest_name = "template"
            help = "template tsv [default: $(DEFAULT_TEMPLATE)]"
    end

    # args.epilog = """
    #     examples:\n
    #     \ua0\ua0 # chloe.jl -t template.tsv -r reference_dir fasta1 fasta2 ...\n
    #     """
    parse_args(ARGS, cmd_args; as_symbols=true)
end

function cmd_main() 
    parsed_args = getargs()
    level = parsed_args[:level]
    Logging.with_logger(Logging.ConsoleLogger(stderr, get(LOGLEVELS, level, Logging.Warn))) do 
        cmd = parsed_args[:_COMMAND_]
        a = parsed_args[cmd]
        if cmd == :gff3
            Sff2Gff.writeallGFF3(;a...)
        elseif cmd == :annotate
            chloe(;a...)
        elseif cmd == :mmap
            SuffixArray.create_mmaps(;a...)
        elseif cmd == :suffix
            SuffixArray.writesuffixarray(;a...)
        end
    end

end

end # module

