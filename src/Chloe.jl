module Chloe

export annotate_batch, annotate, ReferenceDb, AbstractReferenceDb, ChloeConfig, CHLOE_REFS_DIR
export chloe_main


include("annotate_genomes.jl")
include("chloe_cmd.jl")


import .Annotator: annotate_batch, annotate, ReferenceDb, AbstractReferenceDb, ChloeConfig, CHLOE_REFS_DIR
import .CmdLine: chloe_main
end
