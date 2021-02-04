module Chloe

export annotate_one, annotate, readReferences, readDefaultReferences, Reference, Feature
export MayBeIO, MayBeString
export writeallGFF3
export cmd_main
export distributed_main, chloe_distributed, run_broker, get_distributed_args, maybe_launch_broker
export set_global_logger
export annotate_one_task
export readSingleReference, createTargetReference, inverted_repeat

include("ZMQLogger.jl")
include("annotate_genomes.jl")
include("broker.jl")
include("sff2GFF3.jl")
include("WebAPI.jl")
include("chloe_cmd.jl")
include("tasks.jl")

include("chloe_distributed.jl")

import .Sff2Gff: writeallGFF3
# import .ChloeDistributed: distributed_main, chloe_distributed, run_broker, get_distributed_args, maybe_launch_broker
import .Annotator: annotate, annotate_one, readReferences, MayBeIO, MayBeString, Feature
import .CmdLine: cmd_main
import .ZMQLogging: set_global_logger
import .Annotator: readSingleReference, createTargetReference, inverted_repeat
end
