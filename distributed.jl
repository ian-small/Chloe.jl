include("src/ZMQLogger.jl")
include("src/annotate_genomes.jl")
include("src/WebAPI.jl")
include("src/chloe_distributed.jl")
# import .ChloeDistributed: distributed_main

if abspath(PROGRAM_FILE) == @__FILE__
    distributed_main(true)
end
