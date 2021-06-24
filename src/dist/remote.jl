# small thunk for use in
# julia -p 2 -L src/remote.jl
include("../annotate_genomes.jl")
import .Annotator: annotate, annotate_all, annotate_one
