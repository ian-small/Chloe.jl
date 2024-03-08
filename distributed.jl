#!/usr/bin/env julia
if abspath(PROGRAM_FILE) == @__FILE__
    import Chloe
    Chloe.distributed_main(ARGS)
end
