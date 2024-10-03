module Chloe

export annotate, annotate_batch, ReferenceDb, ChloeConfig, AbstractReferenceDb, annotate_one_task
export chloe_main
export distributed_main,  broker_main
export set_global_logger
export ZMQ_ENDPOINT

include("annotate_genomes.jl")

include("dist/ZMQLogger.jl")
include("dist/broker.jl")
include("dist/WebAPI.jl")
include("dist/chloe_cmd.jl")
include("dist/chloe_distributed.jl")

# import .ChloeDistributed: distributed_main, chloe_distributed, run_broker, get_distributed_args, maybe_launch_broker
import .Annotator: annotate_batch, annotate, ReferenceDb, AbstractReferenceDb, ChloeConfig
import .CmdLine: chloe_main
import .Broker: broker_main
import .ZMQLogging: set_global_logger
import .ChloeDistributed: distributed_main, ZMQ_ENDPOINT, annotate_one_task
end
