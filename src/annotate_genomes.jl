module Annotator

using Base: String
using StatsBase: IntegerVector
using XGBoost
export annotate, annotate_one, MayBeIO, MayBeString, ReferenceDb
export read_single_reference!, inverted_repeat, ChloeConfig

import Base

include("utilities.jl")
include("circularity.jl")
include("rotate_genome.jl")
include("hash.jl")
include("align.jl")
include("annotations.jl")
include("orfs.jl")
include("rnas.jl")
include("reference.jl")

import Printf: @sprintf
import JSON
import Crayons: @crayon_str
using StatsBase

const success = crayon"bold green"
const MINIMUM_UNASSIGNED_ORF_LENGTH = Int32(270)

const LACKS_ESSENTIAL_FEATURE = "lacks an essential part of the gene"
const LACKS_START_CODON = "lacks a start codon"
const PREMATURE_STOP_CODON = "has a premature stop codon"
const CDS_NOT_DIVISIBLE_BY_3 = "CDS is not divisible by 3"
const LOW_ANNOTATION_DEPTH = "has low annotation depth, probably spurious"
const INFERIOR_COPY = "probable pseudogene as better copy exists in the genome"
const PARTIAL_IR = "probable pseudogene as overlaps end of inverted repeat"
const OVERLAPPING_FEATURE = "better-scoring feature overlaps with this one"

struct SingleReference
    ref_id::String
    ref_seq::CircularSequence
    ref_features::Union{Nothing,FwdRev{FeatureArray}}
end

datasize(r::SingleReference) = begin
    (sizeof(SingleReference)
     + datasize(r.ref_seq)
     + datasize(r.ref_features)
    )
end

function Base.show(io::IO, r::SingleReference)
    total = datasize(r)
    bp = length(r.ref_seq)
    print(io, "SingleReference[$(r.ref_id)]: $(human(bp))bp, total=$(human(total))")
end

function Base.show(io::IO, r::FwdRev{FwdRev{AlignedBlocks}})
    total = datasize(r)
    ff = length(r.forward.forward)
    fr = length(r.forward.reverse)
    rf = length(r.reverse.forward)
    rr = length(r.reverse.reverse)
    print(io, "Chloë Alignment: [($ff,$fr),($rf,$rr)] total=$(human(total))")
end

MayBeString = Union{Nothing,String}
Strand = Tuple{AAFeature,DFeatureStack}

function flatten(vanno::Vector{Vector{Annotation}})::Vector{Annotation}
    ret = Vector{Annotation}(undef, sum(length(v) for v in vanno; init=0))
    i = 1
    @inbounds for v in vanno
        for a in v
            ret[i] = a
            i += 1
        end
    end
    ret
end

function transfer_annotations(ref_features::FeatureArray, aligned_blocks::BlockTree, src_length::Int32, target_length::Int32)::Vector{Annotation}
    annotations = Vector{Annotation}()
    feature_blocks = intersect(ref_features.feature_tree, aligned_blocks)
    for fb in feature_blocks
        feature = fb[1]
        block = fb[2]
        feature_overlaps = overlaps(range(block.src_index, length=block.blocklength), range(feature.start, length=feature.length), src_length)
        for o in feature_overlaps
            offset5 = o.start - feature.start
            offset3 = o.stop - (feature.start + feature.length - 1)
            if (feature.type == "CDS")
                phase = phase_counter(feature.phase, offset5)
            else
                phase = 0
            end
            # coordinates in target genome
            if o.start > feature.start # alignment must start within feature, so start in target is start of alignment
                tgt_start = block.tgt_index
            else # alignment starts before feature, so start in target is after start of alignment
                tgt_start = mod1(block.tgt_index + circulardistance(block.src_index, feature.start, src_length), target_length)
            end
            push!(annotations, Annotation(ref_features.genome_id, annotation_path(feature), tgt_start, circulardistance(o.start, o.stop, src_length) + 1, offset5, offset3, phase))
        end
    end
    annotations
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

    annotations = flatten(tgt)
    sort!(annotations, by=x -> x.path)
    annotations
end

function score_feature(sff::SFF_Feature, reference_feature_counts::Dict{String,Int}, gmatch::Float32, seq::CircularSequence)
    ref_count = reference_feature_counts[annotation_path(sff.feature)]
    slice = sff.feature.stack[range(sff.feature.start, length=sff.feature.length)]
    sff.stackdepth = Float32((length(slice) > 0 ? sum(slice) : 0) / (ref_count * sff.feature.length))
    sff.relative_length = sff.feature.length / sff.feature.median_length
    sff.gmatch = gmatch
    if sff.feature.type == "CDS"
        codonfrequencies = countcodons(sff.feature, seq)
        sff.coding_prob = xgb_coding_classifier(codonfrequencies)
    end
    sff.feature_prob = feature_xgb(sff.feature.type, sff.feature.median_length, sff.feature.length, sff.stackdepth, gmatch, sff.coding_prob)
end

function fill_feature_stack(target_length::Int32, annotations::Vector{Annotation},
    feature_templates::Dict{String,FeatureTemplate})::Dict{String, CircularVector}
    # assumes annotations are ordered by path
    stacks = Dict{String, CircularVector}()
    for annotation in annotations
        template = get(feature_templates, annotation.path, nothing)
        if template === nothing
            @error "Can't find template for $(annotation.path)"
            continue
        end
        stack = get(stacks, annotation.path, nothing)
        if isnothing(stack)
            stack = stacks[annotation.path] = CircularVector(zeros(Int8, target_length))
        end

        @inbounds for i = annotation.start:annotation.start+annotation.length-one(Int32)
            stack[i] += one(Int8)
        end
    end
    @debug "found ($(length(stacks))) $(human(datasize(stacks))) FeatureStacks from $(length(annotations)) annotations"
    return stacks
end

function align_template(stack::CircularVector, template::FeatureTemplate)::Tuple{Vector{Tuple{Int32,Int}},Int32}
    glen = length(stack)
    median_length = floor(Int32, template.median_length)

    hits = Vector{Tuple{Int32,Int}}()
    score = 0
    @inbounds for nt::Int32 = one(Int32):median_length
        score += stack[nt]
    end

    if score > 0
        push!(hits, (one(Int32), score))
    end

    @inbounds for nt::Int32 = 2:glen
        mt::Int32 = nt + median_length - one(Int32)
        score -= stack[nt-one(Int32)] # remove tail
        score += stack[mt] # add head
        if score > 0
            if isempty(hits)
                push!(hits, (nt, score))
            else
                previous = last(hits)
                if nt ≥ previous[1] + median_length
                    push!(hits, (nt, score))
                elseif score > previous[2]  # if hits overlap, keep the best one
                    pop!(hits)
                    push!(hits, (nt, score))
                end
            end
        end
    end
    # features near start of genome can align twice when template alignment wraps, the check below prevents reporting two hits in this case
    if length(hits) > 1 && last(hits)[1] + median_length - glen > first(hits)[1]
        if last(hits)[2] > first(hits)[2]
            popfirst!(hits)
        else
            pop!(hits)
        end
    end
    isempty(hits) && return hits, median_length
    # sort by descending score
    sort!(hits, by=x -> x[2], rev=true)
    maxscore = hits[1][2]
    # retain all hits scoring over 90% of the the top score
    filter!(x -> (x[2] ≥ maxscore * 0.9 || x[2] ≥ median_length), hits)
    return hits, median_length
end

function do_strand(target_id::String, target_seq::CircularSequence, refs::Vector{SingleReference}, coverages::Dict{String,Float32}, reference_feature_counts::Dict{String,Int},
    strand::Char, blocks_aligned_to_target::Vector{FwdRev{BlockTree}}, feature_templates::Dict{String,FeatureTemplate})::Vector{Vector{SFF_Feature}}

    t4 = time_ns()
    target_length = Int32(length(target_seq))
    annotations = do_annotations(target_id, strand, refs, blocks_aligned_to_target, target_length)

    # strand_feature_stacks is a Dict of feature stacks
    strand_feature_stacks = fill_feature_stack(target_length, annotations, feature_templates)

    t5 = time_ns()
    @info "[$target_id]$strand built feature stacks ($(length(strand_feature_stacks))) $(human(datasize(strand_feature_stacks))): $(ns(t5 - t4))"

    features = Feature[]

    for (path, stack) in strand_feature_stacks
        template = feature_templates[path]
        hits, length = align_template(stack, template)
        isempty(hits) && continue
        for hit in hits
            push!(features, Feature(path, stack, hit[1], length, Int8(0), template.essential, template.median_length))
        end
    end

    gmatch = mean(values(coverages))
    sff_features = Vector{SFF_Feature}(undef, length(features))
    for (i, feature) in enumerate(features)
        refine_boundaries_by_offsets!(feature, annotations, target_length, coverages)
        sff_features[i] = SFF_Feature(feature, 0.0, 0.0, 0.0, 0.0, 0.0)
        score_feature(sff_features[i], reference_feature_counts, gmatch, target_seq)
    end

    # group by feature name on features ordered by mid-point
    target_strand_models::Vector{Vector{SFF_Feature}} = features2models(sort(sff_features, by=x -> x.feature))

    orfs = getallorfs(target_seq, strand, Int32(0))
    # this toys with the feature start, phase etc....
    # refine_gene_models! does not (yet) use the relative_length, stackdepth, feature_prob, coding_prob values but could...
    refine_gene_models!(target_strand_models, target_seq, orfs)

    # this inefficiently updates all features, even those unchanged by refine_gene_models!
    for m in target_strand_models, sf in m
        # update feature data
        score_feature(sf, reference_feature_counts, gmatch, target_seq)
    end

    #=
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
            if orf_frame == f_frame && length(overlaps(range(uorf.start, length = uorf.length), range(f.start, length = f.length), target_length)) > 0
                unassigned = false
            end
        end
        if unassigned
            # count codons
            codonfrequencies = countcodons(uorf, target_seq)
            # predict with XGBCodingClassifier
            coding_prob = xgb_coding_classifier(codonfrequencies)
            push!(target_strand_models, [SFF_Feature(uorf, 0.0, 0.0, gmatch, coding_prob / 2, coding_prob)])
        end
    end
    =#

    @info "[$target_id]$strand built gene models: $(elapsed(t5))"
    return target_strand_models
end

# MayBeIO: write to file (String), IO buffer or create filename based on fasta filename
# Union{IO,String}: read fasta from IO buffer or a file (String)
MayBeIO = Union{String,IO,Nothing}

function align_to_reference(ref::SingleReference, tgt_id::String, tgt_seq::FwdRev{CircularSequence})::Tuple{FwdRev{FwdRev{BlockTree}},Float32}
    start = time_ns()
    ref_length = length(ref.ref_seq)
    tgt_length = length(tgt_seq.forward)

    # align2seqs returns linked list, convert to interval tree
    fchain = align2seqs(ref.ref_seq, tgt_seq.forward, true)
    ff = chain2tree(fchain, ref_length)
    rr = chain2rctree(fchain, ref_length, tgt_length)
    rchain = align2seqs(ref.ref_seq, tgt_seq.reverse, true)
    fr = chain2tree(rchain, ref_length)
    rf = chain2rctree(rchain, ref_length, tgt_length)

    coverage = Float32(100 * target_coverage(fchain, rchain, length(tgt_seq.forward)))

    @info "[$(tgt_id)]± aligned $(ref.ref_id) coverage: $(@sprintf("%.2f", coverage))% $(elapsed(start))"
    # note cross ...
    (FwdRev(FwdRev(ff, rf), FwdRev(fr, rr)), coverage)
end

function inverted_repeat(target::CircularSequence, revtarget::CircularSequence)::AlignedBlock
    fr::AlignedBlocks = ll2vector(align2seqs(target, revtarget, false))
    # sort blocks by length
    fr = sort!(fr, by=b -> b.blocklength, rev=true)
    ir = length(fr) > 0 ? fr[1] : AlignedBlock(0, 0, 0)
    return ir
end

struct ChloeAnnotation
    target_id::String
    target_length::Int32
    coverages::Dict{String,Float32}
    annotation::FwdRev{Vector{SFF_Model}}
end
function annotate_one_worker(db::ReferenceDb,
    target_id::String,
    target::FwdRev{CircularSequence},
    config::ChloeConfig,
)::ChloeAnnotation

    t1 = time_ns()
    target_length = Int32(length(target.forward))

    target_length = length(target.forward)
    @info "[$target_id] seq length: $(target_length)bp"

    refpicks = Vector{Tuple{String,Int}}(undef, 0)

    # find closest gs references
    refhashes = get_gsminhashes(db, config)
    if !isnothing(refhashes)
        hash = minhash(target.forward)
        numrefs = min(config.numgsrefs, length(refhashes))
        append!(refpicks, searchhashes(hash, refhashes)[1:numrefs])
    end

    # find closest chloe references
    refhashes = get_chloeminhashes(db, config)
    if !isnothing(refhashes)
        numrefs = min(config.numchloerefs, length(refhashes))
        append!(refpicks, searchhashes(hash, refhashes)[1:numrefs])
    end
    numrefs = length(refpicks)
    t2 = time_ns()

    @info "[$target_id] picked $numrefs reference(s): $(ns(t2 - t1))"

    blocks_aligned_to_targetf = Vector{FwdRev{BlockTree}}(undef, numrefs)
    blocks_aligned_to_targetr = Vector{FwdRev{BlockTree}}(undef, numrefs)
    coverages = Dict{String,Float32}()
    refs = Vector{SingleReference}(undef, 0)
    reference_feature_counts = Dict{String,Int}()
    # get_single_reference! throws if bad refpick
    refs = [get_single_reference!(db, r[1], reference_feature_counts) for r in refpicks]

    Threads.@threads for i in eachindex(refs)
        a = align_to_reference(refs[i], target_id, target)
        blocks_aligned_to_targetf[i] = a[1].forward
        blocks_aligned_to_targetr[i] = a[1].reverse
        lock(REENTRANT_LOCK) do
            coverages[refpicks[i][1]] = a[2]
        end
    end

    feature_templates = get_templates(db)

    t3 = time_ns()
    @info "[$target_id] aligned: ($(numrefs)) $(human(datasize(blocks_aligned_to_targetf) + datasize(blocks_aligned_to_targetr))) mean coverage: $(geomean(values(coverages))) $(ns(t3 - t2))"

    function watson()
        models = do_strand(target_id, target.forward, refs, coverages, reference_feature_counts,
            '+', blocks_aligned_to_targetf, feature_templates)
        final_models = SFF_Model[]
        for model in filter(m -> !isempty(m), models)
            sff = toSFFModel(feature_templates, model, '+', target.forward, config.sensitivity)
            !isnothing(sff) && push!(final_models, sff)
        end
        final_models
    end

    function crick()
        models = do_strand(target_id, target.reverse, refs, coverages, reference_feature_counts,
            '-', blocks_aligned_to_targetr, feature_templates)
        final_models = SFF_Model[]
        for model in filter(m -> !isempty(m), models)
            sff = toSFFModel(feature_templates, model, '-', target.reverse, config.sensitivity)
            !isnothing(sff) && push!(final_models, sff)
        end
        final_models
    end

    # from https://discourse.julialang.org/t/threads-threads-to-return-results/47382

    sffs_fwd, sffs_rev, ir = fetch.((Threads.@spawn w()) for w in [watson, crick, () -> inverted_repeat(target.forward, target.reverse)])

    ir1 = ir2 = nothing
    if ir.blocklength >= 1000
        @info "[$target_id] inverted repeat $(ir.blocklength)"
        ir1 = SFF_Model("IR-1", 0.0, '+', 1, 1, [SFF_Feature(Feature("IR/repeat_region/1", ir.src_index, ir.blocklength, Int8(0)), 0.0, 0.0, 0.0, 0.0, 0.0)], [])
        ir2 = SFF_Model("IR-2", 0.0, '-', 1, 1, [SFF_Feature(Feature("IR/repeat_region/1", ir.tgt_index, ir.blocklength, Int8(0)), 0.0, 0.0, 0.0, 0.0, 0.0)], [])
    end

    filter_gene_models!(sffs_fwd, sffs_rev, target_length, ir1, ir2)

    if !config.nofilter
        filter!(m -> length(m.warnings) == 0, sffs_fwd)
        filter!(m -> length(m.warnings) == 0, sffs_rev)
    end

    if !isnothing(ir1)
        push!(sffs_fwd, ir1)
        push!(sffs_rev, ir2)
    end
    @info success("[$target_id] Overall: $(elapsed(t1))")
    return ChloeAnnotation(target_id, target_length, coverages, FwdRev(sffs_fwd, sffs_rev))

end

function write_result(result::ChloeAnnotation, asgff3::Bool, output::MayBeIO=nothing)::Tuple{Union{String,IO},String}
    models = update_genecount!(result.annotation)
    ext = asgff3 ? "gff3" : "sff"
    fname = if output !== nothing
        if output isa String
            if isdir(output)
                joinpath(output, "$(result.target_id).$(ext)")
            else
                output # filename
            end
        else
            output # IOBuffer
        end
    else
        "$(result.target_id).$(ext)"
    end
    if !asgff3
        writeSFF(fname, result.target_id, result.target_length, geomean(values(result.coverages)), models)
    else
        writeGFF3(fname, result.target_id, result.target_length, models)
    end

    return fname, result.target_id
end

"""
    annotate_one(refsdir::String, target_id::String, seq::String, [,output_sff_file])

Annotate a single sequence containting a *single* circular
DNA entry

writes an .sff file to `output_sff_file` or uses the sequence id in the
fasta file to write `{seq_id}.sff` in the current directory.

If `output_sff_file` is a *directory* write `{seq_id}.sff` into that
directory.

returns a 2-tuple: (ultimate sff output filename, sequence id)

If `output_sff_file` is an IOBuffer then that buffer will be returned
with the annotation within it.
"""
function annotate_one(db::ReferenceDb,
    target_id::String,
    target::FwdRev{CircularSequence},
    config::ChloeConfig,
    output::MayBeIO=nothing
)::Tuple{Union{String,IO},String}

    result = annotate_one_worker(db, target_id, target, config)
    write_result(result, config.to_gff3, output)

end



function annotate_one(db::ReferenceDb, infile::String, config::ChloeConfig, output::MayBeIO=nothing)
    maybe_gzread(infile) do io
        annotate_one(db, io, config, output)
    end
end

function fasta_reader(infile::IO)::Tuple{String,FwdRev{CircularSequence}}
    reader = FASTA.Reader(infile)
    records = [record for record in reader]
    if isempty(records)
        error("unable to read sequence from $infile")
    elseif length(records) > 2
        @error "$infile contains multiple sequences; Chloë expects a single sequence per file"
    end
    fseq = CircularSequence(FASTA.sequence(LongSequence{DNAAlphabet{4}}, records[1]))
    if length(records) == 2
        rseq = CircularSequence(FASTA.sequence(LongSequence{DNAAlphabet{4}}, records[2]))
    else
        rseq = reverse_complement(fseq)
    end
    if length(fseq) != length(rseq)
        @error "unequal lengths for forward and reverse strands"
    end
    target_id = FASTA.identifier(records[1])
    if isnothing(target_id)
        target_id = "unknown"
    end
    return target_id, FwdRev(fseq, rseq)
end

function annotate_one(db::ReferenceDb, infile::IO, config::ChloeConfig, output::MayBeIO=nothing)
    target_id, seqs = fasta_reader(infile)
    annotate_one(db, target_id, seqs, config, output)
end

function sffname(fafile::String, asgff3::Bool, directory::Union{String,Nothing}=nothing)::String
    ext = asgff3 ? "gff3" : "sff"
    d = if isnothing(directory)
        dirname(fafile)
    else
        directory
    end
    f = basename(fafile)
    base = rsplit(f, '.'; limit=2)[1]

    joinpath(d, "$(base).$(ext)")
end

#= function annotate_one(refsdir::String, infile::String, output::MayBeIO=nothing)

    annotate_one(ReferenceDb(;refsdir=refsdir), infile, ChloeConfig(), output)

end =#

function annotate(db::ReferenceDb, fa_files::Vector{String}, config::ChloeConfig, output::Union{Nothing,String}=nothing)
    n = length(fa_files)
    for infile in fa_files
        maybe_gzread(infile) do io
            annotate_one(db, io, config, if n > 1
                sffname(infile, config.to_gff3, output)
            else
                output
            end)
        end
    end
end

function weighted_mode(values::Vector{Int32}, weights::Vector{Float32})::Int32
    w, v = findmax(StatsBase.addcounts!(Dict{Int32,Float32}(), values, weights))
    return v
end

# uses weighted mode, weighting by alignment length and distance from boundary
function refine_boundaries_by_offsets!(feat::Feature, annotations::Vector{Annotation},
    target_length::Integer, coverages::Dict{String,Float32})
    # grab all the matching features and sort by start in genome of origin to ensure that when iterated in order, most 5' match is found first
    matching_annotations = sort(annotations[findall(x -> x.path == annotation_path(feat), annotations)], by=a -> a.start)
    isempty(matching_annotations) && return #  feat, [], []
    overlapping_annotations = Annotation[]
    minstart = target_length
    maxend = 1

    ## Fix 5' boundary and feature phase
    for annotation in matching_annotations
        annotation.offset5 >= feat.length && continue
        if length(overlaps(range(feat.start, length=feat.length), range(annotation.start, length=annotation.length), target_length)) > 0
            # ignore if we already have an annotation from this genome
            if isnothing(findfirst(a -> a.genome_id == annotation.genome_id, overlapping_annotations))
                minstart = min(minstart, annotation.start)
                push!(overlapping_annotations, annotation)
            end
        end
    end

    n = length(overlapping_annotations)
    end5s = Vector{Int32}(undef, n)
    end5ws = Vector{Float32}(undef, n)
    phases = Vector{Int32}(undef, n)
    @inbounds for (i, annotation) in enumerate(overlapping_annotations)
        # predicted 5' end is annotation start - offset5
        end5s[i] = annotation.start - annotation.offset5
        # weights
        coverage = coverages[annotation.genome_id]
        weight = ((feat.length - (annotation.start - minstart)) / feat.length) * coverage
        end5ws[i] = weight * weight
        phases[i] = Int32(phase_counter(annotation.phase, -annotation.offset5))
    end
    left = feat.start
    if !isempty(end5s)
        left = weighted_mode(end5s, end5ws)
    end
    feat.start = genome_wrap(target_length, left)
    feat.phase = (feat.type ≠ "CDS" || isempty(phases)) ? 0 : weighted_mode(phases, end5ws)

    ## Fix 3' boundary
    empty!(overlapping_annotations)
    for annotation in Iterators.reverse(matching_annotations) # iterate in reverse to ensure 3' annotations are reached first
        annotation.offset3 >= feat.length && continue
        if length(overlaps(range(feat.start, length=feat.length), range(annotation.start, length=annotation.length), target_length)) > 0
            # ignore if we already have an annotation from this genome
            if isnothing(findfirst(a -> a.genome_id == annotation.genome_id, overlapping_annotations))
                maxend = max(maxend, annotation.start + annotation.length - 1)
                push!(overlapping_annotations, annotation)
            end
        end
    end
    n = length(overlapping_annotations)
    end3s = Vector{Int32}(undef, n)
    end3ws = Vector{Float32}(undef, n)
    @inbounds for (i, annotation) in enumerate(overlapping_annotations)
        # predicted 3' end is feature start + feature length + offset3 - 1
        end3s[i] = annotation.start + annotation.length - annotation.offset3 - 1
        # weights
        coverage = coverages[annotation.genome_id]
        weight = ((feat.length - (maxend - (annotation.start + annotation.length - 1))) / feat.length) * coverage
        end3ws[i] = weight * weight
    end
    right = feat.start + feat.length - 1
    if !isempty(end3s)
        right = weighted_mode(end3s, end3ws)
    end
    feat.length = right - left + 1
end

function features2models(sff_features::Vector{SFF_Feature})::Vector{Vector{SFF_Feature}}
    gene_models = [SFF_Feature[]]
    current_model = SFF_Feature[]
    left_border::Int32 = typemax(Int32)
    right_border::Int32 = 0
    for sff_feature in sff_features
        if isempty(current_model)
            push!(current_model, sff_feature)
            left_border = sff_feature.feature.start
            right_border = sff_feature.feature.start + sff_feature.feature.length - 1
            # add feature to model if the gene names match and they are less than 3kb apart; 3kb is an arbitrary limit set to include the longest intron (trnK, ~2.5kb) and may need to be reduced
        elseif current_model[1].feature.gene == sff_feature.feature.gene && max(left_border, sff_feature.feature.start) - min(right_border, sff_feature.feature.start + sff_feature.feature.length - 1) < 3000
            push!(current_model, sff_feature)
            left_border = min(left_border, sff_feature.feature.start)
            right_border = max(right_border, sff_feature.feature.start + sff_feature.feature.length - 1)
        else
            sort!(current_model, by=x -> x.feature.start)
            push!(gene_models, current_model)
            current_model = SFF_Feature[]
            push!(current_model, sff_feature)
            left_border = sff_feature.feature.start
            right_border = sff_feature.feature.start + sff_feature.feature.length - 1
        end
    end
    if length(current_model) > 0
        sort!(current_model, by=x -> x.feature.start)
        push!(gene_models, current_model)
    end
    return gene_models
end

function splice_model(target_seq::CircularSequence, model::Vector{SFF_Feature})::LongDNA{2}
    cds = dna""d    # dynamically allocated so thread-safe
    for sfeat in model
        feat = sfeat.feature
        feat.type == "intron" && continue
        append!(cds, target_seq[feat.start:feat.start+feat.length-1])
    end
    cds = cds[first(model).feature.phase+1:end]
    return cds
end

function startScore(cds::Feature, codon::SubString)
    if codon == "ATG"
        return 1.0
    elseif codon == "GTG"
        if isFeatureName(cds, "rps19")
            return 1.0
        else
            return 0.01
        end
    elseif codon == "ACG"
        if isFeatureName(cds, "ndhD")
            return 1.0
        else
            return 0.01
        end
    else
        return 0
    end
end

function findStartCodon!(cds::Feature, genome_length::Int32, target_seq::CircularSequence)
    # assumes phase has been correctly set
    # search for start codon 5'-3' beginning at cds.start, save result; abort if stop encountered
    start3::Int32 = cds.start + cds.phase
    codon = getcodon(target_seq, start3)
    if isstartcodon(codon, true, true)  # allow ACG and GTG codons if this is predicted start
        cds.start = start3
        cds.length -= cds.phase
        cds.phase = 0
        return cds
    end

    allowACG = false
    allowGTG = false
    if cds.gene == "rps19"
        allowGTG = true
    end

    start3 = 0
    for i::Int32 in cds.start+cds.phase:3:cds.start+cds.phase+genome_length-3
        codon = getcodon(target_seq, i)
        if isstartcodon(codon, allowACG, allowGTG)
            start3 = i >= cds.start + cds.phase ? i : i + genome_length
            break
        end
        if isstopcodon(codon, false)
            break
        end
    end
    # search for start codon 3'-5' beginning at cds.start, save result; abort if stop encountered
    start5::Int32 = 0
    for i::Int32 in cds.start+cds.phase:-3:cds.start+cds.phase-genome_length+3
        codon = getcodon(target_seq, i)
        if isstartcodon(codon, allowACG, allowGTG)
            start5 = i <= cds.start + cds.phase ? i : i - genome_length
            break
        end
        if isstopcodon(codon, false)
            break
        end
    end
    # return cds with start set to nearer of the two choices
    if start3 == 0 && start5 == 0
        # Couldn't find start
        return cds
    end
    if start3 == 0
        cds.length += cds.start - start5
        cds.start = genome_wrap(genome_length, start5)
    elseif start5 == 0
        cds.length += cds.start - start3
        cds.start = genome_wrap(genome_length, start3)
    elseif (start3 - cds.start)^2 < cds.start - start5
        cds.length += cds.start - start3
        cds.start = genome_wrap(genome_length, start3)
    else
        cds.length += cds.start - start5
        cds.start = genome_wrap(genome_length, start5)
    end
    cds.phase = 0
    return cds
end

function setlongestORF!(sfeat::SFF_Feature, orfs::Vector{Feature}, genome_length::Int32)
    # find in-frame and out-of-frame orfs with longest overlap with feature
    # orfs are sorted by descending length
    f = sfeat.feature
    f_frame::Int8 = mod1(f.start + f.phase, 3)
    max_inframe_overlap = 0
    max_outframe_overlap = 0
    inframe_overlap_orf = nothing
    outframe_overlap_orf = nothing
    for orf in orfs
        orf.length < max_inframe_overlap && break
        orf_frame = mod1(orf.start, 3)
        f_frame ≠ orf_frame && orf.length <= max_outframe_overlap && continue
        overlap_array = overlaps(range(orf.start, length=orf.length), range(f.start, length=f.length), genome_length)
        length(overlap_array) == 0 && continue
        overlap = overlap_array[1].stop - overlap_array[1].start + 1
        if f_frame == orf_frame && overlap > max_inframe_overlap
            max_inframe_overlap = overlap
            inframe_overlap_orf = orf
        elseif f_frame ≠ orf_frame && overlap > max_outframe_overlap
            max_outframe_overlap = overlap
            outframe_overlap_orf = orf
        end
    end
    # arbitrary 0.75 proportion threshold for preferring in-frame ORF
    if max_inframe_overlap > max_outframe_overlap || max_inframe_overlap > 0.75 * f.length
        overlap_orf = inframe_overlap_orf
    else
        overlap_orf = outframe_overlap_orf
    end
    if isnothing(overlap_orf)
        f.length = 0
        return
    end
    if overlap_orf.start > f.start # must be an internal stop
        f.phase = 0
        f.length = f.length - (overlap_orf.start - f.start)
        f.start = overlap_orf.start
    else # in case of different frames, have to adjust phase
        f.phase = phase_counter(Int8(0), f.start - overlap_orf.start)
    end
    if sfeat.feature.gene ≠ "rps12A"
        f.length = overlap_orf.start - f.start + overlap_orf.length
    end
end

function refine_boundaries_by_score!(feat1::Feature, feat2::Feature, genome_length::Int32)
    # feat1 should be before feat2
    start = feat1.start + feat1.length - one(Int32)
    range_to_test = min(start, feat2.start):max(start, feat2.start)

    score = sum(feat2.stack[range_to_test]) - sum(feat1.stack[range_to_test])
    maxscore = score
    fulcrum = first(range_to_test) - one(Int32)
    @inbounds for i in range_to_test
        score += 2 * feat1.stack[i]
        score -= 2 * feat2.stack[i]
        if score > maxscore
            maxscore = score
            fulcrum = i
        end
    end
    # fulcrum should point to end of feat1
    feat1.length = fulcrum - feat1.start + one(Int32)
    # fulcrum+1 should point to start of feat2
    feat2.length += feat2.start - (fulcrum + one(Int32))
    feat2.start = genome_wrap(genome_length, fulcrum + one(Int32))
end

function refine_gene_models!(gene_models::Vector{Vector{SFF_Feature}}, target_seq::CircularSequence, orfs::Vector{Feature})

    ## add code to make use of the feature_prob field in the SFF_Features

    genome_length = Int32(length(target_seq))
    last_cds_examined = nothing
    for model in gene_models
        isempty(model) && continue
        # first check feature order and remove out of order features
        pointer = 1
        while true
            pointer == length(model) && break
            sfeature1 = model[pointer]
            sfeature2 = model[pointer+1]
            if sfeature1.feature.order > sfeature2.feature.order # features are out of order
                if sfeature1.feature_prob > sfeature2.feature_prob
                    deleteat!(model, pointer + 1)
                else
                    deleteat!(model, pointer)
                end
                pointer = max(pointer - 1, 1) # back 1 to check order is still OK
            else
                pointer += 1
            end
        end
        ftype = featuretype(model)

        # tRNA or ncRNA, check the ends are complementary
        if ftype == "ncRNA"
            trimRNAends!(target_seq, model)
        end

        if ftype == "tRNA"
            tRNAends!(target_seq, model)
        end

        # if CDS, find phase, stop codon and set feature.length
        if ftype == "CDS"
            local lastexon::SFF_Feature
            for sff in reverse(model)
                if sff.feature.type == "CDS"
                    lastexon = sff
                    break
                end
            end
            setlongestORF!(lastexon, orfs, genome_length)
            last_cds_examined = lastexon.feature
        end

        i = length(model)
        while true
            i == 1 && break
            feature = model[i].feature
            previous_feature = model[i-1].feature
            if feature.length ≤ 0 || (feature.type == "intron" && feature.length < 250)
                deleteat!(model, i)
            elseif previous_feature.length ≤ 0 || (previous_feature.type == "intron" && previous_feature.length < 250)
                deleteat!(model, i - 1)
            else
                gap = feature.start - (previous_feature.start + previous_feature.length)
                if gap ≠ 0 && gap < 100
                    refine_boundaries_by_score!(previous_feature, feature, genome_length)
                    gap = feature.start - (previous_feature.start + previous_feature.length)
                    if feature.length ≤ 0 || (feature.type == "intron" && feature.length < 250)
                        deleteat!(model, i)
                        if i ≤ length(model)
                            feature = model[i].feature
                        end
                    end
                end
                if (gap < 100) && (annotation_path(previous_feature) == annotation_path(feature))
                    previous_feature.length = feature.start - previous_feature.start + feature.length
                    deleteat!(model, i)
                end
            end
            i -= 1
        end

        isempty(model) && continue
        # if CDS, find start codon and set feature.start
        first_exon = first(model).feature
        if first_exon.type == "CDS"
            if lastexon.feature.gene ≠ "rps12B"
                findStartCodon!(first_exon, genome_length, target_seq)
            end
        end
    end
end

function get_model_boundaries(model::SFF_Model, glength::Int32)::UnitRange{Int32}
    genestart = first(model.features).feature.start
    geneend = last(model.features).feature.start + last(model.features).feature.length - 1
    length::Int32 = geneend > genestart ? geneend - genestart + 1 : glength + geneend - genestart + 1
    return range(genestart, length=length)
end

function mean_stackdepth(model::SFF_Model)::Float64
    sum = 0.0
    for sff in model.features
        sum += sff.stackdepth
    end
    return sum / length(model.features)
end

function toSFFModel(feature_templates::Dict{String,FeatureTemplate}, model::Vector{SFF_Feature}, strand::Char, target_seq::CircularSequence, sensitivity::Real)::Union{Nothing,SFF_Model}
    gene = first(model).feature.gene
    type = featuretype(model)
    type == "" && return nothing
    warnings = String[]
    expected_exons = String[]
    for n in 1:5
        possible_exon = gene * "/" * type * "/" * string(n)
        if haskey(feature_templates, possible_exon)
            push!(expected_exons, possible_exon)
        end
    end
    exon_count = 0
    model_prob = 0.0
    for sff in model
        if sff.feature.type ≠ "intron"
            model_prob += sff.feature_prob
            exon_count += 1
        end
    end
    if exon_count < length(expected_exons)
        push!(warnings, LACKS_ESSENTIAL_FEATURE)
    end

    if startswith(gene, "unassigned_orf")
        model_prob = first(model).feature_prob
    else
        model_prob = model_prob / length(expected_exons)
    end
    exceeds_sensitivity = false
    if model_prob ≥ sensitivity || isnan(model_prob)
        exceeds_sensitivity = true
    end

    if exceeds_sensitivity && (type == "CDS")
        cds = splice_model(target_seq, model)
        if gene ≠ "rps12B" && !isstartcodon(getcodon(cds, Int32(1)), true, true)
            push!(warnings, LACKS_START_CODON)
        end
        for i::Int32 in 1:3:length(cds)-2
            if isstopcodon(getcodon(cds, i))
                push!(warnings, PREMATURE_STOP_CODON)
                break
            end
        end
        if length(cds) % 3 ≠ 0
            push!(warnings, CDS_NOT_DIVISIBLE_BY_3)
        end
    end
    return exceeds_sensitivity ? SFF_Model(gene, model_prob, strand, 1, exon_count, model, warnings) : nothing
end

function write_model2SFF(outfile::IO, model::SFF_Model)
    isnothing(model) && return
    model_id = model.gene * "/" * string(model.gene_count)
    for sff in model.features
        f = sff.feature
        id = "$(model.gene)/$(model.gene_count)/$(f.type)/$(f.order)"
        write(outfile, id)
        write(outfile, "\t")
        write(outfile, join([model.strand, string(f.start), string(f.length), string(f.phase)], "\t"))
        write(outfile, "\t")
        write(outfile, join([@sprintf("%.3g", sff.relative_length), @sprintf("%.3g", sff.stackdepth), @sprintf("%.3g", sff.gmatch),
                @sprintf("%.3g", sff.feature_prob), @sprintf("%.3g", sff.coding_prob)], "\t"))
        write(outfile, "\t")
        write(outfile, join(model.warnings, "; "))
        write(outfile, "\n")
    end
    for warning in model.warnings
        @warn("$(model_id) $(warning)")
    end
end

function calc_maxlengths(models::FwdRev{Vector{Vector{SFF_Model}}})::Dict{String,Int32}
    maxlengths = Dict{String,Int32}()
    function add_model(models)
        for model in models
            isempty(model) && continue
            m = first(model)
            maxlengths[m.gene] = max(m.gene_length, get(maxlengths, m.gene, 0))
        end
    end

    add_model(models.forward)
    add_model(models.reverse)
    maxlengths
end

coding_xgb_model = Booster(model_file=joinpath(@__DIR__, "coding_xgb.model"))
noncoding_xgb_model = Booster(model_file=joinpath(@__DIR__, "noncoding_xgb.model"))
const MAXFEATURELENGTH = 7000
function feature_xgb(ftype::String, median_length::Float32, featurelength::Int32, fdepth::Float32, gmatch::Float32, codingprob::Float32)::Float32
    featurelength ≤ 0 && return Float32(0.0)
    fdepth ≤ 0 && return Float32(0.0)
    rtlength = median_length / MAXFEATURELENGTH
    rflength = featurelength / MAXFEATURELENGTH
    frtlength = featurelength / median_length
    match = gmatch / 100
    local pred
    if ftype == "CDS"
        pred = XGBoost.predict(coding_xgb_model, [rtlength rflength frtlength fdepth match codingprob])
    else
        pred = XGBoost.predict(noncoding_xgb_model, [rtlength rflength fdepth match])
    end
    return pred[1]
end

# only some gene_models are allowed to overlap
function allowed_model_overlap(m1, m2)::Bool
    if m1.gene == "trnK-UUU" && m2.gene == "matK"
        return true
    end
    if m1.gene == "matK" && m2.gene == "trnK-UUU"
        return true
    end
    if m1.gene == "ndhC" && m2.gene == "ndhK"
        return true
    end
    if m1.gene == "ndhK" && m2.gene == "ndhC"
        return true
    end
    if m1.gene == "psbC" && m2.gene == "psbD"
        return true
    end
    if m1.gene == "psbD" && m2.gene == "psbC"
        return true
    end
    if m1.gene == "atpB" && m2.gene == "atpE"
        return true
    end
    if m1.gene == "atpE" && m2.gene == "atpB"
        return true
    end
    if m1.gene == "rps3" && m2.gene == "rpl22"
        return true
    end
    if m1.gene == "rpl22" && m2.gene == "rps3"
        return true
    end
    return false
end

function filter_gene_models!(fwd_models::Vector{SFF_Model}, rev_models::Vector{SFF_Model}, glength::Int32, ir1::Union{Nothing,SFF_Model}, ir2::Union{Nothing,SFF_Model})

    # hard floor on stackdepth
    for model in fwd_models
        if mean_stackdepth(model) < 0.01
            push!(model.warnings, LOW_ANNOTATION_DEPTH)
        end
    end
    for model in rev_models
        if mean_stackdepth(model) < 0.01
            push!(model.warnings, LOW_ANNOTATION_DEPTH)
        end
    end

    # filter out models of genes where better model of same gene exists
    # filter out models of genes where better model of overlapping gene exists
    function warningcheck!(models)
        for (i, model1) in enumerate(models)
            for j in i+1:length(models)
                model2 = models[j]
                bmodel1 = get_model_boundaries(model1, glength)
                bmodel2 = get_model_boundaries(model2, glength)
                if model1.gene == model2.gene
                    warning = INFERIOR_COPY
                elseif length(overlaps(bmodel1, bmodel2, glength)) > 0 && !allowed_model_overlap(model1, model2)
                    warning = OVERLAPPING_FEATURE
                else
                    continue
                end
                if model1.gene_prob > model2.gene_prob ##no tolerance for unequal probs
                    push!(model2.warnings, warning)
                elseif model2.gene_prob > model1.gene_prob
                    push!(model1.warnings, warning)
                elseif mean_stackdepth(model1) > mean_stackdepth(model2) # if probs are equal, decide on stackdepth
                    push!(model2.warnings, warning)
                elseif mean_stackdepth(model2) > mean_stackdepth(model1)
                    push!(model1.warnings, warning)
                end
            end
        end
    end
    warningcheck!(fwd_models)
    warningcheck!(rev_models)

    #deal with duplicated modes on opposite strands, e.g. in the IR
    for model1 in fwd_models, model2 in rev_models
        if model1.gene == model2.gene
            model1_boundaries = get_model_boundaries(model1, glength)
            model2_boundaries = get_model_boundaries(model2, glength)
            if !isnothing(ir1) && !isnothing(ir2)
                #if both models are within the IRs, both can be kept, if not, one is flagged
                ir1_boundaries = get_model_boundaries(ir1, glength)
                ir2_boundaries = get_model_boundaries(ir2, glength)
                if model1.strand == '+'
                    model1_maxIRintersect = max(length(circularintersect(model1_boundaries, ir1_boundaries, glength)), length(circularintersect(model1_boundaries, reverse_complement(ir2_boundaries, glength), glength)))
                else
                    model1_maxIRintersect = max(length(circularintersect(model1_boundaries, reverse_complement(ir1_boundaries, glength), glength)), length(circularintersect(model1_boundaries, ir2_boundaries, glength)))
                end
                if model2.strand == '+'
                    model2_maxIRintersect = max(length(circularintersect(model2_boundaries, ir1_boundaries, glength)), length(circularintersect(model2_boundaries, reverse_complement(ir2_boundaries, glength), glength)))
                else
                    model2_maxIRintersect = max(length(circularintersect(model2_boundaries, reverse_complement(ir1_boundaries, glength), glength)), length(circularintersect(model2_boundaries, ir2_boundaries, glength)))
                end
                if (model1_maxIRintersect == 0 || model2_maxIRintersect == 0) # one or other model outside the IRs, deal with as INFERIOR_COPY
                    if model1.gene_prob > model2.gene_prob ##no tolerance for unequal probs
                        push!(model2.warnings, INFERIOR_COPY)
                    elseif model2.gene_prob > model1.gene_prob
                        push!(model1.warnings, INFERIOR_COPY)
                    elseif mean_stackdepth(model1) > mean_stackdepth(model2) # if probs are equal, decide on stackdepth
                        push!(model2.warnings, INFERIOR_COPY)
                    elseif mean_stackdepth(model2) > mean_stackdepth(model1)
                        push!(model1.warnings, INFERIOR_COPY)
                    end
                elseif model1_maxIRintersect == length(model1_boundaries) && model2_maxIRintersect == length(model2_boundaries)  # both models inside the IRs, no warnings
                    continue
                else   # overlap IR, keep best model
                    # find IR intersect of feature in model1 that overlaps end of IR
                    local f::Feature, frange::UnitRange{Int32}, intersect::UnitRange{Int32}
                    for sff in model1.features
                        f = sff.feature
                        frange = range(f.start, length = f.length)
                        if model1.strand == '+'
                            intersect = circularintersect(frange, ir1_boundaries, glength)
                            if length(intersect) == 0
                                intersect = circularintersect(frange, reverse_complement(ir2_boundaries, glength), glength)
                            end
                        else
                            intersect = circularintersect(frange, reverse_complement(ir1_boundaries, glength), glength)
                            if length(intersect) == 0
                                intersect = circularintersect(frange, ir2_boundaries, glength)
                            end
                        end
                        0 < length(intersect) < f.length && break
                    end
                    # find relative depth in non-IR region of feature stack
                    rd1 = sum(f.stack[setdiff(frange, intersect)])
                    # find feature in model2 that overlaps end of IR
                    for sff in model2.features
                        f = sff.feature
                        frange = range(f.start, length = f.length)
                        if model2.strand == '+'
                            intersect = circularintersect(frange, ir1_boundaries, glength)
                            if length(intersect) == 0
                                intersect = circularintersect(frange, reverse_complement(ir2_boundaries, glength), glength)
                            end
                        else
                            intersect = circularintersect(frange, reverse_complement(ir1_boundaries, glength), glength)
                            if length(intersect) == 0
                                intersect = circularintersect(frange, ir2_boundaries, glength)
                            end
                        end
                        0 < length(intersect) < f.length && break
                    end
                    # find relative depth in non-IR region of featute stack
                    rd2 = sum(f.stack[setdiff(frange, intersect)])
                    if rd1 > rd2
                        push!(model2.warnings, PARTIAL_IR)
                    elseif rd2 > rd1
                        push!(model1.warnings, PARTIAL_IR)
                    end
                end
            end
        end
    end
end

function modelID!(model_ids::Dict{String,Int32}, model::SFF_Model)
    gene_name = model.gene
    instance_count::Int8 = 1
    model_id = "$(gene_name)_$(instance_count)"
    while get(model_ids, model_id, 0) ≠ 0
        instance_count += 1
        model_id = "$(gene_name)_$(instance_count)"
    end
    model_ids[model_id] = instance_count
    model.gene_count = instance_count
end

function update_genecount!(models::FwdRev{Vector{SFF_Model}})::FwdRev{Vector{SFF_Model}}
    model_ids = Dict{String,Int32}()
    fwd = SFF_Model[]
    rev = SFF_Model[]
    for model in models.forward
        isnothing(model) && continue
        isempty(model.features) && continue
        modelID!(model_ids, model) # updates model.gene_count
        push!(fwd, model)

    end
    for model in models.reverse
        isnothing(model) && continue
        isempty(model.features) && continue
        modelID!(model_ids, model)# updates model.gene_count
        push!(rev, model)

    end
    return FwdRev(fwd, rev)
end

function writeSFF(outfile::Union{String,IO},
    id::String, # NCBI id
    genome_length::Int32,
    mean_coverage::Float32,
    models::FwdRev{Vector{SFF_Model}})

    function out(outfile::IO)
        write(outfile, id, "\t", string(genome_length), "\t", @sprintf("%.3f", mean_coverage), "\n")
        for model in models.forward
            write_model2SFF(outfile, model)
        end
        for model in models.reverse
            write_model2SFF(outfile, model)
        end
    end
    if outfile isa String
        maybe_gzwrite(outfile::String) do io
            out(io)
        end
    else
        out(outfile)
    end
end

function sff2gffcoords(f::Feature, strand::Char, genome_length::Integer)
    start = f.start
    finish = f.start + f.length - 1
    length = finish - start + 1
    finish = mod1(finish, genome_length)

    if strand == '-'
        start = genome_length - finish + 1
        finish = start + length - 1
        start = mod1(start, genome_length)
        finish = mod1(finish, genome_length)
    end
    return (start, finish, length)
end

function writeGFF3(outfile::Union{String,IO},
    genome_id::String, # NCBI id
    genome_length::Int32,
    models::FwdRev{Vector{SFF_Model}})

    function header(outfile::IO)
        write(outfile, "##gff-version 3.2.1\n")
        write(outfile, join([genome_id, "Chloe", "region", '1', string(genome_length), ".", "+", "0", "Is_circular=true"], "\t"))
        write(outfile, "\n")
        write(outfile, join([genome_id, "Chloe", "source", '1', string(genome_length), ".", "+", "1", "ID=source-null;gbkey=source"], "\t"))
        write(outfile, "\n")
        write(outfile, "###\n")
    end

    function out(outfile::IO)
        header(outfile)
        allmodels = sort(vcat(models.forward, models.reverse), by=m -> sff2gffcoords(first(m.features).feature, m.strand, genome_length)[1])
        for model in allmodels
            write_model2GFF3(outfile, model, genome_id, genome_length)
        end
    end

    if outfile isa String
        maybe_gzwrite(outfile::String) do io
            out(io)
        end
    else
        out(outfile)
    end
end

function merge_adjacent_features!(model::SFF_Model)
    f1_index = 1
    f2_index = 2
    while f2_index <= length(model.features)
        f1 = model.features[f1_index].feature
        f2 = model.features[f2_index].feature
        # if adjacent features are same type, merge them into a single feature
        if f1.type == f2.type && f2.start - f1.start - f1.length ≤ 100
            @debug "[$(genome_id)]$(strand) merging adjacent $(f1.path) and $(f2.path)"
            f1.length = f2.start - f1.start + f2.length
            deleteat!(model.features, f2_index)
        else
            f1_index += 1
            f2_index += 1
        end
    end
end

function featuretype(model::Vector{SFF_Feature})
    type = ""
    for sff in model
        if sff.feature.type ≠ "intron"
            type = sff.feature.type
            break
        end
    end
    return type
end

function write_model2GFF3(outfile, model::SFF_Model, genome_id::String, genome_length::Int32)

    function write_line(type, start, finish, pvalue, id, parent, phase="."; key="Parent")
        l = [genome_id, "Chloe", type, start, finish, @sprintf("%.3e", pvalue), model.strand, phase]
        write(outfile, join(l, "\t"))
        write(outfile, "\t", "ID=", id, ";", key, "=", parent, "\n")
    end

    merge_adjacent_features!(model)

    id = model.gene
    if model.gene_count > 1
        id = id * "-" * string(model.gene_count)
    end

    start = minimum(f.feature.start for f in model.features)
    ft = featuretype(model.features)
    if ft == "CDS" && !startswith(id, "rps12A")
        last(model.features).feature.length += 3 # add stop codon
    end
    finish = maximum(f.feature.start + f.feature.length - 1 for f in model.features)
    length = finish - start + 1

    if model.strand == '-'
        start = genome_length - finish + 1
        if start < 1
            start += genome_length
        end
        finish = start + length - 1
    end

    start = mod1(start, genome_length)
    finish = mod1(finish, genome_length)

    # gene
    if startswith(model.gene, "IR")
        write_line("repeat_region", start, finish, model.gene_prob, id, model.gene; key="Name")
        write(outfile, "###\n")
        return
    else
        write_line("gene", start, finish, 1.0 - model.gene_prob, id, model.gene; key="Name")
    end
    # RNA product
    parent = id
    #= if ft == "CDS"
        parent =  id * ".mRNA"
        write_line("mRNA", start, finish, 1.0 - model.gene_prob, parent, id)
    elseif ft == "rRNA"
        parent = id * ".rRNA"
        write_line("rRNA", start, finish, 1.0 - model.gene_prob, parent, id)
    elseif ft == "tRNA"
        parent =  id * ".tRNA"
        write_line("tRNA", start, finish, 1.0 - model.gene_prob, parent, id)
    end =#
    # exons
    for sff in model.features
        f = sff.feature
        type = f.type
        #= if type == "tRNA" || type == "rRNA"
            type = "exon"
        end =#
        start, finish, length = sff2gffcoords(f, model.strand, genome_length)

        phase = type == "CDS" ? string(f.phase) : "."
        write_line(type, start, finish, 1.0 - sff.feature_prob, id * "." * type * "." * string(f.order), parent, phase)
    end
    write(outfile, "###\n")
end

end # module
