#!/usr/bin/env julia
if abspath(PROGRAM_FILE) == @__FILE__
    include("src/dist/ZMQLogger.jl")
    include("src/annotate_genomes.jl")
    include("src/dist/WebAPI.jl")
    include("src/dist/broker.jl")
    include("src/dist/chloe_distributed.jl")
    # import .ChloeDistributed: distributed_main
    distributed_main(true)
end
