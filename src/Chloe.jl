module Chloe

export annotate_batch, annotate, ReferenceDb, AbstractReferenceDb, ChloeConfig
export chloe_main


include("annotate_genomes.jl")
include("chloe_cmd.jl")


import .Annotator: annotate_batch, annotate, ReferenceDb, AbstractReferenceDb, ChloeConfig
import .CmdLine: chloe_main
end
