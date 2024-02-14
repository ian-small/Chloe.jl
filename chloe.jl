include("src/annotate_genomes.jl")
include("src/dist/chloe_cmd.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    import .CmdLine: cmd_main
    cmd_main()
end
