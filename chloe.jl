include("src/annotate_genomes.jl")
include("src/sff2GFF3.jl")
include("src/SuffixArray.jl")
include("src/chloe_cmd.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    import .CmdLine: cmd_main
    cmd_main()
end
