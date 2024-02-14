#!/usr/bin/env julia
project = joinpath(dirname(@__FILE__), "..", "Project.toml")
import Pkg; 
# Pkg.activate(project)
(open(project) |> Pkg.TOML.parse)["deps"] |> keys |> collect |> Pkg.add
