#!/usr/bin/env julia
# use this if you have installed Chloe as a package.
# Upside: can make use of julia's package pre-compilation.
# Downside: you can't add new workers to a running server with this
# setup (use distributed.jl). YMMV
using Chloe
using Distributed
args = get_distributed_args()
args = maybe_launch_broker(args)
addprocs(args[:workers], topology=:master_worker)
@everywhere workers() using Chloe
chloe_distributed(false;args...)