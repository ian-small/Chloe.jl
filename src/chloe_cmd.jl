
module CmdLine
export cmd_main

import ArgParse: ArgParseSettings, @add_arg_table!, parse_args

import ..Annotator
import ..Sff2Gff
import ..SuffixArray

include("globals.jl")
include("rotate_genome.jl")


function chloe(;refsdir=DEFAULT_REFS, fasta_files=String[],
    template=DEFAULT_TEMPLATE, output::Union{Nothing,String}=nothing,
    forward_only::Bool=false, verbose::Bool=true)
    if refsdir == "default"
        refsdir = normpath(joinpath(HERE, "..", DEFAULT_REFS))
    end
    if template == "default"
        template = normpath(joinpath(HERE, "..", DEFAULT_TEMPLATE))
    end
    Annotator.annotate(refsdir, template, fasta_files, output; forward_only=forward_only,
        verbose=verbose)

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
        "rotate"
            help = "rotate circular genomes to standard form"
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
        "--forward-only", "-o"
            action = :store_true
            help = "only save forward"
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
        "--forward-only"
            action = :store_true
            help = "only use forward sequences"
    end

    @add_arg_table! cmd_args["rotate"] begin
        "fasta-files"
            arg_type = String
            nargs = '+'
            required = true
            action = :store_arg
            help = "fasta file to process"
        "--flip-SSC", "-S"
            action = :store_true
            help = "flip orientation of small single-copy region"
        "--flip-LSC", "-L"
            action = :store_true
            help = "flip orientation of large single-copy region"
        "--extend", "-e"
            default = 0
            arg_type = Int
            help = "add n bases from start to end of sequence to allow mapping to wrap ends [use -1 for maximum extent]"
        "--directory", "-d"
            arg_type = String
            help = "output directory to place rotated fasta files (use '-' to write to stdout)" 
        "--width", "-w"
            arg_type = Int
            default = 80
            help = "line width of fasta output"
    end

    # args.epilog = """
    #     examples:\n
    #     \ua0\ua0 # chloe.jl -t template.tsv -r reference_dir fasta1 fasta2 ...\n
    #     """
    parse_args(ARGS, cmd_args; as_symbols=true)
end

function cmd_main() 
    parsed_args = getargs()
    level = lowercase(parsed_args[:level])
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
        elseif cmd == :rotate
            rotateGenome(;a...)
        end
    end

end

end # module

