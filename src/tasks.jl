
# put these in the global namespace
import .ZMQLogging: annotation_local_storage, set_global_logger, TASK_KEY
import .Annotator: annotate_one, MayBeIO, MayBeString, readReferences, Reference


#### these are only used by chloe_distributed ####

function annotate_one_task(fasta::MayBeString, output::MayBeIO, task_id::MayBeString)
    annotation_local_storage(TASK_KEY, task_id)
    try
        # the global REFERENCE should have been
        # sent to the worker process by main process
        @debug "using $(Main.REFERENCE)"
        annotate_one(Main.REFERENCE, fasta, output)
    finally
        annotation_local_storage(TASK_KEY, nothing)
    end
end


function annotate_one_task(fasta::Union{String,IO}, task_id::MayBeString)
    annotation_local_storage(TASK_KEY, task_id)
    try
        # the global REFERENCE should have been
        # sent to the worker process by main process
        @debug "using $(Main.REFERENCE)"
        annotate_one(Main.REFERENCE, fasta, IOBuffer())
    finally
        annotation_local_storage(TASK_KEY, nothing)
    end
end
