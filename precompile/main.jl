# include("annotate_genomes.jl")
module Chloe
using ArgParse
const DEBUG = get(ENV, "DEBUG", "false") == "true"
 
function main(;refsdir = "reference_1116", fasta_files = String[], verbose = false,
    template = "optimised_templates.v2.tsv")
    println("refs: " * refsdir)
    for f in fasta_files
        println("fasta: " * f)
    end
end


# const julia_v07 = VERSION > v"0.7-"

s = ArgParseSettings(autofix_names = true)  # turn "-" into "_" for arg names.

@add_arg_table! s begin
    "fasta-files"
        arg_type = String
        nargs = '+'
        required = true
        action = :store_arg
        help = "fasta files to process"
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

end
s.epilog = """
    examples:\n
    \ua0\ua0 # main.jl -t template.tsv -r reference_dir fasta1 fasta2 ...\n
    """

function real_main() 
    parsed_args = parse_args(ARGS, s; as_symbols = true)
    # filter!((k, v)->v ∉ (nothing, false), parsed_args)
    filter!(kv->kv.second ∉ (nothing, false), parsed_args)
    println(parsed_args)
    main(;parsed_args...)
end

Base.@ccallable function julia_main()::Cint
    try
        real_main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    real_main()
end

end # module