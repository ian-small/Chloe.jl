module Chloe

export annotate_one, annotate, Feature
export MayBeIO, MayBeString
export writeallGFF3
export cmd_main
export distributed_main, chloe_distributed, run_broker, get_distributed_args, maybe_launch_broker
export set_global_logger
export annotate_one_task
export read_single_reference!, inverted_repeat

include("dist/ZMQLogger.jl")
include("annotate_genomes.jl")
include("dist/broker.jl")
include("dist/WebAPI.jl")
include("dist/chloe_cmd.jl")
include("dist/tasks.jl")

include("dist/chloe_distributed.jl")

# import .ChloeDistributed: distributed_main, chloe_distributed, run_broker, get_distributed_args, maybe_launch_broker
import .Annotator: annotate, annotate_one, MayBeIO, MayBeString, Feature, ReferenceDb
import .CmdLine: cmd_main
import .ZMQLogging: set_global_logger
import .Annotator: read_single_reference!, inverted_repeat
end
