
# put these in the global namespace
import .ZMQLogging: annotation_local_storage, set_global_logger, TASK_KEY
import .Annotator: annotate_one, MayBeIO, MayBeString, CholeConfig


#### these are only used by chloe_distributed ####

function annotate_one_task(fasta::String, output::MayBeIO, task_id::MayBeString, config::ChloeConfig)
    annotation_local_storage(TASK_KEY, task_id)
    try
        annotate_one(Main.REFERENCE, fasta, config, output)
    finally
        annotation_local_storage(TASK_KEY, nothing)
    end
end


function annotate_one_task(fasta::Union{String,IO}, task_id::MayBeString, config::ChloeConfig)
    annotation_local_storage(TASK_KEY, task_id)
    try
        annotate_one(Main.REFERENCE, fasta, config, IOBuffer())
    finally
        annotation_local_storage(TASK_KEY, nothing)
    end
end
