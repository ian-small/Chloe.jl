#!/usr/bin/env julia
if abspath(PROGRAM_FILE) == @__FILE__
    include("src/ZMQLogger.jl")
    include("src/annotate_genomes.jl")
    include("src/WebAPI.jl")
    include("src/broker.jl")
    include("src/chloe_distributed.jl")
    # import .ChloeDistributed: distributed_main
    distributed_main(true)
end
