# small thunk for use in
# julia -p 2 -L src/dist/remote.jl
using Pkg
Pkg.activate(pwd())
include("../annotate_genomes.jl")
import .Annotator: annotate, annotate_one
