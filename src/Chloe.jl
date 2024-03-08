module Chloe

export annotate, annotate_batch, ReferenceDbFromDir, ReferenceDb, ChloeConfig, AbstractReferenceDb
# export MayBeIO, MayBeString
export chloe_main
export distributed_main, chloe_distributed, run_broker, broker_main, get_distributed_args, maybe_launch_broker
export set_global_logger
export annotate_one_task
export ZMQ_CLIENT

include("annotate_genomes.jl")

include("dist/ZMQLogger.jl")
include("dist/broker.jl")
include("dist/WebAPI.jl")
include("dist/chloe_cmd.jl")
include("dist/tasks.jl")
include("dist/chloe_distributed.jl")

# import .ChloeDistributed: distributed_main, chloe_distributed, run_broker, get_distributed_args, maybe_launch_broker
import .Annotator: annotate_batch, annotate, ReferenceDb, ReferenceDbFromDir, AbstractReferenceDb, ChloeConfig
import .CmdLine: chloe_main
import .Broker: broker_main
import .ZMQLogging: set_global_logger
end
