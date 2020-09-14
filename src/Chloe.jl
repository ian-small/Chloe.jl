module Chloe

export start_broker, annotate_one, annotate, readReferences
export create_mmaps, writesuffixarray
export writeallGFF3

include("annotate_genomes.jl")
include("broker.jl")
include("sff2GFF3.jl")
include("make_suffix_array.jl")

end
