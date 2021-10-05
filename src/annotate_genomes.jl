module Annotator

export annotate, annotate_one, read_references, Reference, MayBeIO, MayBeString
export read_single_reference, inverted_repeat

import Base

include("Utilities.jl")
include("rotate_genome.jl")
include("hash.jl")
include("AlignGenomes.jl")
include("Annotations.jl")
include("orfs.jl")

# import Dates: Time, Nanosecond

import Printf:@sprintf
import JSON
import Crayons:@crayon_str
using StatsBase
using GLM

const success = crayon"bold green"
const MINIMUM_UNASSIGNED_ORF_LENGTH = Int32(270)

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

function read_single_reference!(refdir::String, refID::AbstractString, reference_feature_counts::Dict{String, Int})::Union{SingleReference, Nothing}
    try
        ref = FASTA.Record()
        if !isdir(refdir); refdir = dirname(refdir); end
        path = findfastafile(refdir, refID)
        if isnothing(path)
            msg = "unable to find $(refID) fasta file in $(refdir)!"
            @error msg
            throw(ArgumentError(msg))
        end
        reader = open(FASTA.Reader, path)
        read!(reader, ref)
        close(reader)
        path = normpath(joinpath(refdir, refID * ".sff"))
        if !isfile(path)
            msg = "unable to find $(path)!"
            @error msg
            throw(ArgumentError(msg))
        end
        ref_features = read_features!(path, reference_feature_counts)
        return SingleReference(refID, CircularSequence(FASTA.sequence(ref)), ref_features)
    catch e
        return nothing
    end
end

MayBeString = Union{Nothing,String}
Strand = Tuple{AAFeature,DFeatureStack}

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

function do_annotations(target_id::String, strand::Char, refs::Vector{SingleReference}, blocks_aligned_to_target::Vector{FwdRev{BlockTree}}, target_length::Int32)

    function do_one(refsrc, ref_features, blocks)
        st = time_ns()
        annotations = transfer_annotations(ref_features.forward, blocks.forward, length(refsrc), target_length)
        annotations = vcat(annotations, transfer_annotations(ref_features.reverse, blocks.reverse, length(refsrc), target_length))
        @debug "[$(target_id)]$(strand) $(refsrc)± overlaps $(length(annotations)): $(elapsed(st))"
        return annotations
    end
    tgt = Vector{Vector{Annotation}}(undef, length(refs))

    Threads.@threads for i in eachindex(refs)
        blocks = blocks_aligned_to_target[i]
        tgt[i] = do_one(refs[i].ref_seq, refs[i].ref_features, blocks)
    end

    # annotations = collect(Iterators.flatten(tgt))
    annotations = flatten(tgt)
    sort!(annotations, by=x -> x.path)
    annotations
end

function score_feature(sff::SFF_Feature, maxtemplatelength::Float32, stack::FeatureStack, reference_feature_counts::Dict{String,Int}, gmatch::Float32, seq::CircularSequence)
    ref_count = reference_feature_counts[annotation_path(sff.feature)]
    sff.stackdepth = Float32(sum(stack.stack[range(sff.feature.start, length=sff.feature.length)]) / (ref_count * sff.feature.length))
    sff.relative_length = sff.feature.length / stack.template.median_length
    sff.gmatch = gmatch
    sff.feature_prob = feature_xgb(maxtemplatelength, stack.template, sff.relative_length, sff.stackdepth, gmatch)
    if sff.feature.type == "CDS"
        sff.coding_prob = glm_coding_classifier(countcodons(sff.feature, seq))
    end
end

function do_strand(target_id::String, target_seq::CircularSequence, refs::Vector{SingleReference}, coverages::Dict{String,Float32}, reference_feature_counts::Dict{String,Int},
    strand::Char, blocks_aligned_to_target::Vector{FwdRev{BlockTree}}, feature_templates::Dict{String,FeatureTemplate}, maxtemplatelength::Float32)::Vector{Vector{SFF_Feature}}

    # println(strand, '\t', blocks_aligned_to_target)

    t4 = time_ns()
    target_length = Int32(length(target_seq))
    annotations = do_annotations(target_id, strand, refs, blocks_aligned_to_target, target_length)

    #println(strand, '\t', annotations)

    # strand_feature_stacks is basically grouped by annotations.path
    strand_feature_stacks = fill_feature_stack(target_length, annotations, feature_templates)

    t5 = time_ns()
    @info "[$target_id]$strand built feature stacks ($(length(strand_feature_stacks))) $(human(datasize(strand_feature_stacks))): $(ns(t5 - t4))"

    # orfmap = threeframes(target_seq)
    features = Feature[]
    path_to_stack = Dict{String,FeatureStack}()

    for stack in strand_feature_stacks
        hits, length = align_template(stack)
        isempty(hits) && continue
        for hit in hits
            push!(features, Feature(stack.path, hit[1], length, 0))
        end
        if haskey(path_to_stack, stack.path)
            @error "duplicate feature path: $(stack.path)"
        end
        path_to_stack[stack.path] = stack
    end

    #println(strand, '\t', features)

    gmatch = mean(values(coverages))
    sff_features = Vector{SFF_Feature}(undef, length(features))
    for (i, feature) in enumerate(features)
        refine_boundaries_by_offsets!(feature, annotations, target_length, coverages)
        sff_features[i] = SFF_Feature(feature, 0.0, 0.0, 0.0, 0.0, 0.0)
        stack = path_to_stack[annotation_path(feature)]
        score_feature(sff_features[i], maxtemplatelength, stack, reference_feature_counts, gmatch, target_seq)        
    end

    #println(strand, '\t', sff_features)

    # group by feature name on features ordered by mid-point
    target_strand_models::Vector{Vector{SFF_Feature}} = features2models(sort(sff_features, by = x -> x.feature))

    #println(strand, '\t', target_strand_models)

    orfs = getallorfs(target_seq, strand, Int32(0))
    # this toys with the feature start, phase etc....
    # refine_gene_models! does not (yet) use the relative_length, stackdepth, feature_prob, coding_prob values but could...
    refine_gene_models!(target_strand_models, target_seq, path_to_stack, orfs)

    # this inefficiently updates all features, even those unchanged by refine_gene_models!
    for m in target_strand_models, sf in m
        # update feature data
        f = sf.feature
        stack = path_to_stack[annotation_path(f)]
        score_feature(sf, maxtemplatelength, stack, reference_feature_counts, gmatch, target_seq) 
    end

    #println(strand, '\t', target_strand_models)

    # add any unassigned orfs
    for uorf in orfs
        uorf.length < MINIMUM_UNASSIGNED_ORF_LENGTH && continue
        # println(uorf)
        unassigned = true
        orf_frame = mod1(uorf.start, 3)
        for m in target_strand_models, sf in m
            f = sf.feature
            # println(f)
            f_frame = mod1(f.start + f.phase, 3)
            if orf_frame == f_frame && length(overlaps(range(uorf.start, length=uorf.length), range(f.start, length=f.length), target_length)) > 0
                unassigned = false
            end
        end
        if unassigned
            # count codons
            codonfrequencies = countcodons(uorf, target_seq)
            # predict with GLMCodingClassifier
            coding_prob = glm_coding_classifier(codonfrequencies)
            push!(target_strand_models, [SFF_Feature(uorf, 0.0, 0.0, gmatch, coding_prob / 2, coding_prob)])
        end
    end

    @info "[$target_id]$strand built gene models: $(elapsed(t5))"
    return target_strand_models
end

# MayBeIO: write to file (String), IO buffer or create filename based on fasta filename
# Union{IO,String}: read fasta from IO buffer or a file (String)
MayBeIO = Union{String,IO,Nothing}

function align_to_reference(ref::SingleReference, tgt_id::String, tgt_seq::CircularSequence, rev_tgt_seq::CircularSequence, mask::Union{Nothing,FwdRev{CircularMask}}=nothing)::Tuple{FwdRev{FwdRev{BlockTree}},Float32}
    start = time_ns()
    ref_length = length(ref.ref_seq)
    tgt_length = length(tgt_seq)

    # align2seqs returns linked list, convert to interval tree
    fchain = align2seqs(ref.ref_seq, tgt_seq, mask.forward)
    ff = chain2tree(fchain, ref_length)
    rr = chain2rctree(fchain, ref_length, tgt_length)
    rchain = align2seqs(ref.ref_seq, rev_tgt_seq, mask.reverse)
    fr = chain2tree(rchain, ref_length)
    rf = chain2rctree(rchain, ref_length, tgt_length)

    coverage = Float32(100 * target_coverage(fchain, rchain, length(tgt_seq)))

    @info "[$(tgt_id)]± aligned $(ref.ref_id) coverage: $(@sprintf("%.2f", coverage))% $(elapsed(start))"
    # note cross ...
    (FwdRev(FwdRev(ff, rf), FwdRev(fr, rr)), coverage)
end

function inverted_repeat(target::CircularSequence, revtarget::CircularSequence)::AlignedBlock
    fr::AlignedBlocks = ll2vector(align2seqs(target, revtarget;))
    # sort blocks by length
    fr = sort!(fr, by=b -> b.blocklength, rev=true)
    ir = length(fr) > 0 ? fr[1] : AlignedBlock(0, 0, 0)
    return ir
end

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

`reference` are the reference annotations (see `read_references`)
"""
function annotate_one(refsdir::String, numrefs::Int, refhashes::Union{Nothing,Dict{String,Vector{Int64}}}, feature_templates::Dict{String,FeatureTemplate}, sensitivity::Float16, target_id::String,
     target_forward_strand::CircularSequence, target_reverse_strand::CircularSequence, output::MayBeIO, gff3::Bool, nofilter::Bool)::Tuple{Union{String,IO},String}

    t1 = time_ns()

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

    # sanity checks
    target_length = length(target_forward_strand)
    n = count(isambiguous, target_forward_strand.sequence)
    r = n / target_length
    if r > .01
        @warn("sequence [$(target_id)] contains too many ambiguous nucleotides: $(@sprintf "%.1f" r * 100)% ")
        return fname, target_id
    end
    n = count(isgap, target_forward_strand.sequence)
    if n > 0
        @warn("sequence [$(target_id)] contains gaps which will be removed before analysis")
        ungap!(target_forward_strand.sequence)
        ungap!(target_reverse_strand.sequence)
        target_length = Int32(length(target_forward_strand))
    end
    
    @info "[$target_id] seq length: $(target_length)bp"

    # find best references
    emask = entropy_mask(CircularSequence(target_forward_strand.sequence), Int32(KMERSIZE))
    if !isnothing(refhashes)
        nmaskedtarget = copy(target_forward_strand.sequence)
        for (i, bit) in enumerate(emask)
            if bit; nmaskedtarget[i] = DNA_N; end
        end
        hash = minhash(nmaskedtarget, KMERSIZE, SKETCHSIZE)
        numrefs = min(numrefs, length(refhashes))
        refpicks = searchhashes(hash, refhashes)[1:numrefs]
    else
        numrefs = 1
        refpicks = [(split(basename(refsdir), '.')[1], 0)]
    end

    t2 = time_ns()

    @info "[$target_id] picked $(numrefs) reference(s): $(ns(t2 - t1))"

    blocks_aligned_to_targetf = Vector{FwdRev{BlockTree}}(undef, numrefs)
    blocks_aligned_to_targetr = Vector{FwdRev{BlockTree}}(undef, numrefs)
    coverages = Dict{String,Float32}()
    refs = Vector{SingleReference}(undef, 0)
    masks = FwdRev(emask, CircularMask(reverse(emask.m)))
    reference_feature_counts = Dict{String, Int}()

    picks_to_remove = []
    for i in 1:numrefs
        ref = read_single_reference!(refsdir, refpicks[i][1], reference_feature_counts)
        if isnothing(ref)
            push!(picks_to_remove, i)
        else
            push!(refs, ref)
        end
    end
    deleteat!(refpicks, picks_to_remove)
    
    Threads.@threads for i in eachindex(refs)
        a = align_to_reference(refs[i], target_id, target_forward_strand, target_reverse_strand, masks)
        blocks_aligned_to_targetf[i] = a[1].forward
        blocks_aligned_to_targetr[i] = a[1].reverse
        lock(REENTRANT_LOCK)
        coverages[refpicks[i][1]] = a[2]
        unlock(REENTRANT_LOCK)
    end

    t3 = time_ns()
    @info "[$target_id] aligned: ($(numrefs)) $(human(datasize(blocks_aligned_to_targetf) + datasize(blocks_aligned_to_targetr))) mean coverage: $(geomean(values(coverages))) $(ns(t3 - t2))" 

    maxtemplatelength = maximum([t.median_length for t in values(feature_templates)])

    function watson()
        models = do_strand(target_id, target_forward_strand, refs, coverages, reference_feature_counts,
            '+', blocks_aligned_to_targetf, feature_templates, maxtemplatelength)
        final_models = SFF_Model[]
        for model in filter(m -> !isempty(m), models)
            sff = toSFF(feature_templates, model, '+', target_forward_strand, sensitivity)
            !isnothing(sff) && push!(final_models, sff)
        end
        final_models
    end
        
    function crick()
        models = do_strand(target_id, target_reverse_strand, refs, coverages, reference_feature_counts,
            '-', blocks_aligned_to_targetr, feature_templates, maxtemplatelength)
        final_models = SFF_Model[]
        for model in filter(m -> !isempty(m), models)
            sff = toSFF(feature_templates, model, '-', target_reverse_strand, sensitivity)
            !isnothing(sff) && push!(final_models, sff)
        end
        final_models
    end

    # from https://discourse.julialang.org/t/threads-threads-to-return-results/47382

    sffs_fwd, sffs_rev, ir = fetch.((Threads.@spawn w()) for w in [watson, crick, () -> inverted_repeat(target_forward_strand, target_reverse_strand)]) 

    filter_gene_models!(sffs_fwd, sffs_rev, target_length)
    if !nofilter
        filter!(m -> length(m.warnings) == 0, sffs_fwd)
        filter!(m -> length(m.warnings) == 0, sffs_rev)
    end
   
    if ir.blocklength >= 1000
        @info "[$target_id] inverted repeat $(ir.blocklength)"
        push!(sffs_fwd, SFF_Model("IR-1", 0.0, '+', 1, 1, [SFF_Feature(Feature("IR/repeat_region/1", ir.src_index, ir.blocklength, 0), 0.0, 0.0, 0.0, 0.0, 0.0)], []))
        push!(sffs_rev, SFF_Model("IR-2", 0.0, '-', 1, 1, [SFF_Feature(Feature("IR/repeat_region/1", ir.tgt_index, ir.blocklength, 0), 0.0, 0.0, 0.0, 0.0, 0.0)], []))
    else
        ir = nothing
    end
    
    writeSFF(fname, target_id, target_length, geomean(values(coverages)), FwdRev(sffs_fwd, sffs_rev))
    if gff3
        writeGFF3(join(split(fname, ".")[1:end - 1], ".") * ".gff3", target_id, target_length, FwdRev(sffs_fwd, sffs_rev))
    end
    
    @info success("[$target_id] Overall: $(elapsed(t1))")
    return fname, target_id

end
"""
    annotate_one(reference::Reference, fasta::Union{String,IO})

Annotate a fasta file. Maybe a file name or an IOBuffer
returns a 2-tuple (sff annotation as an IOBuffer, sequence id)

`reference` are the reference annotations (see `read_references`)
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

function annotate(refsdir::String, numrefs::Int, hashfile::Union{String,Nothing}, templates_file::String, fa_files::Vector{String}, sensitivity::Float16, output::MayBeString, gff3::Bool, nofilter::Bool)
    
    numrefs == 1 ? refhashes = nothing : refhashes = readminhashes(hashfile)

    feature_templates = read_templates(templates_file)

    for infile in fa_files
        reader = open(FASTA.Reader, infile)
        records = Vector{FASTA.Record}()
        for record in reader
            push!(records, record)
        end
        if isempty(records)
            @error "unable to read sequence from $infile"
        elseif length(records) > 2
            @error "$infile contains multiple sequences; Chloë expects a single sequence per file"
        end
        fseq = CircularSequence(FASTA.sequence(records[1]))
        if length(records) == 2
            rseq = CircularSequence(FASTA.sequence(records[2]))
        else
            rseq = reverse_complement(fseq)
        end
        if length(fseq) != length(rseq)
            @error "unequal lengths for forward and reverse strands"
        end
        annotate_one(refsdir, numrefs, refhashes, feature_templates, sensitivity, FASTA.identifier(records[1]), fseq, rseq, output, gff3, nofilter)
        close(reader)
    end
end

end # module
