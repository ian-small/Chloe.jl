module Annotator

export annotate, annotate_one, readReferences, Reference, MayBeIO, MayBeString
export readSingleReference, createTargetReference, inverted_repeat

import Base

include("Utilities.jl")
include("hash.jl")
include("AlignGenomes.jl")
include("Annotations.jl")


# import Dates: Time, Nanosecond

import Printf: @sprintf
import JSON
import Crayons: @crayon_str

const success = crayon"bold green"
const NUM_REFS = 16

struct SingleReference
    ref_id::String
    ref_seq::CircularSequence
    ref_features::Union{Nothing,FwdRev{FeatureArray}}
end

datasize(r::SingleReference) = begin
    (sizeof(Reference)
    + datasize(r.refsrc)
    + datasize(r.ref_features)
    )
end

function Base.show(io::IO, r::SingleReference)
    total = datasize(r)
    bp = datasize(length(r.refsrc))
    print(io, "SingleReference[$(r.refsrc)]: $(human(bp))bp, total=$(human(total))")
end

function Base.show(io::IO, r::FwdRev{FwdRev{AlignedBlocks}})
    total = datasize(r)
    ff = length(r.forward.forward)
    fr = length(r.forward.reverse)
    rf = length(r.reverse.forward)
    rr = length(r.reverse.reverse)
    print(io, "Chloë Alignment: [($ff,$fr),($rf,$rr)] total=$(human(total))")
end

function verify_refs(refsdir, template)
    # used by master process to check reference directory
    # *before* starting worker processes...

    # TODO: read json file and really check...
    if !isdir(refsdir)
        msg = "Reference directory $(refsdir) is not a directory!"
        @error msg
        throw(ArgumentError(msg))
    end

    for f in [template, joinpath(refsdir, "ReferenceOrganisms.json")]
        if !isfile(f)
            msg = "missing file: $f"
            @error msg
            throw(ArgumentError(msg))
        end
    end
    files = readdir(refsdir)    
    sff = findall(x -> endswith(x, ".sff"), files)
    if length(sff) == 0
        msg = "no reference .sff files in $(refsdir)!"
        @error msg
        throw(ArgumentError(msg))
    end
end

function readSingleReference(refdir::String, refID::String)::SingleReference
    ref = FASTA.Record()
    path = findfastafile(refdir, refID)
    if isnothing(path)
        msg = "unable to find $(refID) fasta file in $(refdir)!"
        @error msg
        throw(ArgumentError(msg))
    end
    reader = open(FASTA.Reader, path)
    read!(reader, ref)
    close(reader)
    ref_features = readFeatures(normpath(joinpath(refdir, refID * ".sff")))
    SingleReference(refID, CircularSequence(FASTA.sequence(ref)), ref_features)
end

MayBeString = Union{Nothing,String}
Strand = Tuple{AAFeature,DFeatureStack}
AAlignedBlocks = Vector{FwdRev{AlignedBlocks}}

function flatten(vanno::Vector{Vector{Annotation}})::Vector{Annotation}
    ret = Vector{Annotation}(undef, sum(length(v) for v in vanno))
    i = 1
    @inbounds for v in vanno
        for a in v
            ret[i] = a
            i += 1
        end
    end
    ret
end

function do_annotations(target_id::String, strand::Char, refs::Vector{SingleReference}, blocks_aligned_to_target::AAlignedBlocks)
    # this takes about 4secs!
    function do_one(refsrc, ref_features, blocks)
        st = time_ns()
        annotations = findOverlaps(ref_features.forward, blocks.forward)
        annotations = vcat(annotations, findOverlaps(ref_features.reverse, blocks.reverse))
        @debug "[$(target_id)]$(strand) $(refsrc)± overlaps $(length(annotations)): $(elapsed(st))"
        return annotations
    end
    tgt = Vector{Vector{Annotation}}(undef, NUM_REFS)

    #Threads.@threads for i in 1:NUM_REFS
    for i in 1:NUM_REFS
        blocks = blocks_aligned_to_target[i]
        tgt[i] = do_one(refs[i].ref_seq, refs[i].ref_features, blocks)
    end

    # annotations = collect(Iterators.flatten(tgt))
    annotations = flatten(tgt)
    sort!(annotations, by = x -> x.path)
    annotations
end

function do_strand(target_id::String, target_seq::CircularSequence, start_ns::UInt64, refs::Vector{SingleReference}, coverages::Dict{String,Float32},
    strand::Char, blocks_aligned_to_target::AAlignedBlocks, feature_templates::Dict{String,FeatureTemplate})::Strand

    annotations = do_annotations(target_id, strand, refs, blocks_aligned_to_target)

    t4 = time_ns()
    @info "[$target_id]$strand overlapping ref annotations ($(length(annotations))) $(human(datasize(annotations))): $(ns(t4 - start_ns))"

    # strand_feature_stacks is basically grouped by annotations.path
    target_length = Int32(length(target_seq))
    strand_feature_stacks, shadow = fillFeatureStack(target_length, annotations, feature_templates)

    t5 = time_ns()
    @info "[$target_id]$strand ref features stacks ($(length(strand_feature_stacks))) $(human(datasize(strand_feature_stacks))): $(ns(t5 - t4))"

    features = Vector{Feature}()
    # sizehint!(fa, length(strand_feature_stacks))
    path_to_stack = Dict{String,FeatureStack}()
    # so... features will be ordered by annotations.path
    for stack in strand_feature_stacks
        left_border, length = alignTemplateToStack(stack, shadow)
        left_border == 0 && continue
        depth, coverage = getDepthAndCoverage(stack, left_border, length)
        if ((depth >= stack.template.threshold_counts) && (coverage >= stack.template.threshold_coverage))
            push!(features, Feature(stack.path, left_border, length, 0))
            if haskey(path_to_stack, stack.path)
                @error "duplicate feature path: $(stack.path)"
            end
            path_to_stack[stack.path] = stack
        else
            @debug "[$target_id]$strand below threshold: $(stack.path): depth=$(@sprintf "%.3f" depth) coverage=$(@sprintf "%.3f" coverage)"
        end
    end

    t6 = time_ns()
    @info "[$target_id]$strand aligning templates ($(length(features))): $(ns(t6 - t5))"

    for feature in features
        refineMatchBoundariesByOffsets!(feature, annotations, target_length, coverages)
    end

    t7 = time_ns()
    @info "[$target_id]$strand refining match boundaries: $(ns(t7 - t6))"

    # group by feature name on **ordered** features getFeatureName()
    target_strand_models = groupFeaturesIntoGeneModels(features)
    # this toys with the feature start, phase etc....
    target_strand_models = refineGeneModels!(target_strand_models, target_seq, annotations, path_to_stack)

    @info "[$target_id]$strand refining gene models: $(elapsed(t7))"
    return target_strand_models, path_to_stack
end

# MayBeIO: write to file (String), IO buffer or create filename based on fasta filename
# Union{IO,String}: read fasta from IO buffer or a file (String)
MayBeIO = Union{String,IO,Nothing}

function createTargetReference(fasta::Union{String,IO})::SingleReference
    if fasta isa String && endswith(fasta, ".mmap")
        refSAs, refRAs, refloops = read_mmap_suffix(fasta)
        target_length = length(refSAs.forward)
        _, f = splitpath(fasta)
        refsrc = splitext(f)[1]
        return SingleReference(refsrc,
            refloops,
            refSAs,
            refRAs,
            nothing, # no .sff file
            target_length,
            true, false)
    end
    seq_id, seq = readFasta(fasta)
    createTargetReference(seq_id, seq)
end

function createTargetReference(target_id::String, target_seqf::LongDNASeq)::SingleReference
    target_seqr = revComp(target_seqf)
    target_length = Int32(length(target_seqf))
    
    targetloopf = target_seqf * target_seqf[1:end - 1]
    targetloopr = target_seqr * target_seqr[1:end - 1]

    targetloopf = MappedPtrString(targetloopf)
    targetloopr = MappedPtrString(targetloopr)
    
    target_saf = makeSuffixArray(targetloopf, true)
    target_raf = makeSuffixArrayRanksArray(target_saf)

    target_sar = makeSuffixArray(targetloopr, true)
    target_rar = makeSuffixArrayRanksArray(target_sar)

    SingleReference(target_id,
        FwdRev(targetloopf, targetloopr),
        FwdRev(target_saf, target_sar),
        FwdRev(target_raf, target_rar),
        nothing, # no .sff file
        target_length,
        false, false)
end

function align(ref::SingleReference, tgt_id::String, tgt_seq::CircularSequence, rev_tgt_seq::CircularSequence)::Tuple{FwdRev{FwdRev{AlignedBlocks}},Float32}
    start = time_ns()

    #expected = Dict("JX512022" => 188, "NC_002202" => 515, "Z00044" => 504, "AP000423" => 574, "NC_030504" => 508, "NC_005086" => 509, "NC_024542" => 446, "NC_001666" => 573,
    #                "NC_031333" => 485, "NC_026040" => 302, "KT634228" => 335, "NC_016986" => 313, "NC_004543" => 347, "NC_001319" => 594, "MF177093" => 163, "AP005672" => 546)

    #align2seqs returns linked list
    #but subsequent functions expecting Vector
    ff::AlignedBlocks = ll2vector(align2seqs(ref.ref_seq.sequence, tgt_seq.sequence))
    rr = revCompBlocks(ff,length(ref.ref_seq.sequence), length(tgt_seq.sequence))
    
    fr::AlignedBlocks = ll2vector(align2seqs(ref.ref_seq.sequence,rev_tgt_seq.sequence))
    rf = revCompBlocks(fr,length(ref.ref_seq.sequence), length(tgt_seq.sequence))
    
    # if expected[ref.ref_id] ≠ length(ff)
    # #if ref.ref_id == "NC_001666"
    #      println(ref.ref_id,ff)
    # end

    coverage = Float32(100 * target_coverage(ff, rf, length(tgt_seq)))

    @info "[$(tgt_id)]± aligned $(ref.ref_id) coverage: $(@sprintf("%.2f", coverage))% $(elapsed(start))"
    # note cross ...
    (FwdRev(FwdRev(ff, rf), FwdRev(fr, rr)), coverage)
end

function inverted_repeat(target::CircularSequence, revtarget::CircularSequence)::AlignedBlock
    fr::AlignedBlocks = ll2vector(align2seqs(target.sequence,revtarget.sequence))
    # sort blocks by length
    fr = sort!(fr, by= b -> b.blocklength, rev=true)
    ir = length(fr) > 0 ? fr[1] : AlignedBlock(0, 0, 0)
    return ir
end

function avg_coverage(target::SingleReference, a::FwdRev{AlignedBlocks})
    coverage = blockCoverage(a.forward) + blockCoverage(a.reverse)
    coverage /= (target.target_length * 2)
    coverage
end

function avg_coverage(target::SingleReference, a::FwdRev{FwdRev{AlignedBlocks}})
    avg_coverage(target, a.forward)
end

Models = Vector{Vector{SFF}}

"""
    annotate_one(references::Reference, seq_id::String, seq::String, [,output_sff_file])

Annotate a single sequence containting a *single* circular
DNA entry

writes an .sff file to `output_sff_file` or uses the sequence id in the
fasta file to write `{seq_id}.sff` in the current directory.

If `output_sff_file` is a *directory* write `{seq_id}.sff` into that
directory.

returns a 2-tuple: (ultimate sff output filename, sequence id)

If `output_sff_file` is an IOBuffer then that buffer will be returned
with the annotation within it

`reference` are the reference annotations (see `readReferences`)
"""
function annotate_one(refsdir::String, refhashes::Dict{String,Vector{Int64}}, templates_file::String, target_id::String, target_seq::CircularSequence, output::MayBeIO)::Tuple{Union{String,IO},String}

    t1 = time_ns()

    # sanity checks
    target_length = length(target_seq)
    n = count(isambiguous, target_seq.sequence)
    r = n / target_length
    if r > .1
        error("sequence [$(target_id)] contains too many ambiguous nucleotides: $(@sprintf "%.1f" r * 100)% ")
    end
    n = count(isgap, target_seq.sequence)
    if n > 0
        error("sequence [$(target_id)] contains gaps which will be removed before analysis")
        target_seq = ungap!(target_seq)
        target_length = Int32(length(target_seq))
    end

    rev_target_seq = reverse_complement(target_seq)
    
    @info "[$target_id] seq length: $(target_length)bp"

    # find best references
    hash = minhash(target_seq.sequence, KMERSIZE, SKETCHSIZE)
    refpicks = searchhashes(hash, refhashes)[1:NUM_REFS]

    t2 = time_ns()

    @info "[$target_id] picked $(NUM_REFS) references: $(ns(t2 - t1))"

    blocks_aligned_to_targetf = AAlignedBlocks(undef, NUM_REFS)
    blocks_aligned_to_targetr = AAlignedBlocks(undef, NUM_REFS)
    coverages = Dict{String,Float32}()
    refs = Vector{SingleReference}(undef, NUM_REFS)

    Threads.@threads for i in 1:NUM_REFS
    #for i in 1:NUM_REFS
        ref = readSingleReference(refsdir, refpicks[i][1])
        a = align(ref, target_id, target_seq, rev_target_seq)
        blocks_aligned_to_targetf[i] = a[1].forward
        blocks_aligned_to_targetr[i] = a[1].reverse
        lock(REENTRANT_LOCK)
        coverages[refpicks[i][1]] = a[2]
        unlock(REENTRANT_LOCK)
        refs[i] = ref
    end

    t3 = time_ns()
    
    @info "[$target_id] aligned: ($(NUM_REFS)) $(human(datasize(blocks_aligned_to_targetf) + datasize(blocks_aligned_to_targetr))) $(ns(t3 - t2))" 

    feature_templates = readTemplates(templates_file)

    function watson()
        models, stacks = do_strand(target_id, target_seq, t3, refs, coverages,
            '+', blocks_aligned_to_targetf, feature_templates)

        [toSFF(model, target_seq, stacks) 
            for model in filter(m -> !isempty(m), models)]
    end

    function crick()
        models, stacks = do_strand(target_id, rev_target_seq, t3, refs, coverages,
        '-', blocks_aligned_to_targetr, feature_templates)

        [toSFF(model, rev_target_seq, stacks) 
            for model in filter(m -> !isempty(m), models)]
    end

    # from https://discourse.julialang.org/t/threads-threads-to-return-results/47382

    sffs_fwd, sffs_rev, ir = fetch.((Threads.@spawn w()) for w in [watson, crick, () -> inverted_repeat(target_seq, rev_target_seq)]) 
    
    # sffs_fwd = watson()
    # sffs_rev = crick()
    # ir = inverted_repeat(target_seq, rev_target_seq)
   
    if ir.blocklength >= 1000
        @info "[$target_id] inverted repeat $(ir.blocklength)"
    else
        ir = nothing
    end
    
    fname = if output !== nothing
        if output isa String
            if isdir(output)
                joinpath(output, "$(target_id).sff")
            else
                output # filename
            end
        else
            output # IOBuffer
        end
    else
        "$(target_id).sff"
    end
    writeSFF(fname, target_id, target_length, FwdRev(sffs_fwd, sffs_rev), ir)

    @info success("[$target_id] Overall: $(elapsed(t1))")
    return fname, target_id

end
"""
    annotate_one(reference::Reference, fasta::Union{String,IO})

Annotate a fasta file. Maybe a file name or an IOBuffer
returns a 2-tuple (sff annotation as an IOBuffer, sequence id)

`reference` are the reference annotations (see `readReferences`)
"""
# function annotate_one(reference::Reference, fasta::Union{String,IO})::Tuple{Union{String,IO},String}
#     target_id, target_seqf = readFasta(fasta)
#     annotate_one(reference, target_id, target_seqf, IOBuffer())
# end
# function annotate_one(reference::Reference, fasta::Union{String,IO}, output::MayBeIO)::Tuple{Union{String,IO},String}
#     target_id, target_seqf = readFasta(fasta)
#     annotate_one(reference, target_id, target_seqf, output)
# end

# function annotate_all(reference::Reference, fasta::Union{String,IO})
#     for (target_id, target_seqf) in iterFasta(fasta)
#         annotate_one(reference, target_id, target_seqf, nothing)
#     end
# end

function annotate(refsdir::String, hashfile::String, templates::String, fa_files::Vector{String}, output::MayBeString;
    verbose::Bool=true)

    refhashes = readminhashes(hashfile)

    for infile in fa_files
        reader = open(FASTA.Reader, infile)
        record = FASTA.Record()
        while !eof(reader)
            read!(reader, record)
            annotate_one(refsdir, refhashes, templates, FASTA.identifier(record), CircularSequence(FASTA.sequence(record)), output)
        end
        close(reader)
    end
end

end # module
