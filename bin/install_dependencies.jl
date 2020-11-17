#!/usr/bin/env julia
using Pkg
#println("Activating environment in $(pwd())")
pkg"activate ."
println("Installing packages..."); flush(stdout);
pkg"instantiate"
pkg"precompile"
println("Done!")
