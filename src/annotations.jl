import Printf: @sprintf
import StatsBase
using BioSequences
using XGBoost

mutable struct Feature <: AbstractInterval{Int32}
    gene::String
    type::String
    order::UInt8
    start::Int32
    length::Int32
    # phase is the number of nucleotides to skip at the start of the sequence 
    # to be in the correct reading frame
    phase::Int8
    function Feature(path, s, l, p)
        this = new()
        tags = split(path, "/")
        this.gene = tags[1]
        this.type = length(tags) == 4 ? tags[3] : tags[2]
        this.order = length(tags) == 4 ? parse(UInt8, tags[4][1]) : parse(UInt8, tags[3][1])
        this.start = s
        this.length = l
        this.phase = p
        this
    end
end

Base.first(f::Feature) = f.start
function Base.last(f::Feature)
    f.start + f.length - one(Int32)
end

function Base.isless(a::Feature, b::Feature)
    return a.gene == b.gene ? isless(a.start + a.length / 2, b.start + b.length / 2) : isless(a.gene, b.gene)
end

datasize(f::Feature) = sizeof(Feature) + sizeof(f.annotation_path()) + sum(sizeof(p) for p in f._path_components)

const annotation_path(f::Feature) = begin
    join([f.gene, f.type, string(f.order)], "/")
end

AFeature = Vector{Feature}
AAFeature = Vector{AFeature}

struct FeatureTree
    ftree::IntervalTree{Int32,Feature}
    length::Int32
    wrapped_intervals::Dict{Feature,Feature}
end

function Base.push!(ctree::FeatureTree, interval::Feature)
    if last(interval) ≤ ctree.length
        push!(ctree.ftree, interval)
    else
        first_interval = Feature(annotation_path(interval), first(interval), ctree.length - first(interval) + one(Int32), zero(Int8))
        second_interval = Feature(annotation_path(interval), 1, last(interval) - ctree.length, zero(Int8))
        push!(ctree.ftree, first_interval)
        ctree.wrapped_intervals[first_interval] = interval
        push!(ctree.ftree, second_interval)
        ctree.wrapped_intervals[second_interval] = interval
    end
end

function Base.intersect(ftree::FeatureTree, btree::BlockTree)::Vector{Pair{Feature,AlignedBlock}}
    intervals = collect(intersect(ftree.ftree, btree.btree))
    return unique(map(x -> Pair(get(ftree.wrapped_intervals, x[1], x[1]), get(btree.wrapped_intervals, x[2], x[2])), intervals))
end

# entire set of Features for one strand of one genome
struct FeatureArray
    genome_id::String
    genome_length::Int32
    strand::Char
    feature_tree::FeatureTree
end

# extended Feature struct to add prediction info
mutable struct SFF_Feature
    feature::Feature
    relative_length::Float32
    stackdepth::Float32
    gmatch::Float32
    feature_prob::Float32
    coding_prob::Float32
end

function read_sff_features!(file::String, reference_feature_counts::Dict{String,Int})::FwdRev{FeatureArray}
    open(file) do f
        header = split(readline(f), '\t')
        genome_id = header[1]
        genome_length = parse(Int32, header[2])
        r_features = FeatureTree(IntervalTree{Int32,Feature}(), genome_length, Dict{Feature,Feature}())
        f_features = FeatureTree(IntervalTree{Int32,Feature}(), genome_length, Dict{Feature,Feature}())
        while !eof(f)
            fields = split(readline(f), '\t')
            startswith(fields[1], "IR") && continue # don't use IR annotations
            startswith(fields[1], "unassigned") && continue # don't use unassigned annotations
            length(fields) ≥ 11 && length(fields[11]) > 0 && continue # don't use annotations with warnings
            feature = Feature(fields[1], parse(Int, fields[3]), parse(Int, fields[4]), parse(Int, fields[5]))
            path = annotation_path(feature)
            count = get(reference_feature_counts, path, nothing)
            if isnothing(count)
                reference_feature_counts[path] = 1
            else
                reference_feature_counts[path] = count + 1
            end
            if fields[2][1] == '+'
                push!(f_features, feature)
            else
                push!(r_features, feature)
            end
        end
        f_strand_features = FeatureArray(genome_id, genome_length, '+', f_features)
        r_strand_features = FeatureArray(genome_id, genome_length, '-', r_features)
        return FwdRev(f_strand_features, r_strand_features)
    end
end

# part or all of a Feature annotated by alignment
struct Annotation
    genome_id::String
    path::String
    start::Int32
    length::Int32
    # offsets are the distance from the edge of the annotation to the end of 
    # the original feature
    offset5::Int32
    offset3::Int32
    # phase is the number of nucleotides to skip at the start of the sequence 
    # to be in the correct reading frame
    phase::Int8
end
datasize(a::Annotation) = sizeof(Annotation) + sizeof(a.path)# genome_id is shared + sizeof(a.genome_id)
datasize(v::Vector{Annotation}) = sum(datasize(a) for a in v)

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

struct FeatureTemplate
    path::String  # similar to .sff path
    essential::Bool
    median_length::Float32 # median length of feature
end

datasize(s::FeatureTemplate) = sizeof(FeatureTemplate) + sizeof(s.path)

function read_templates(file::String)::Dict{String,FeatureTemplate}
    if !isfile(file)
        error("\"$(file)\" is not a file")
    end
    if filesize(file) === 0
        error("no data in \"$(file)!\"")
    end
    templates = Dict{String,FeatureTemplate}()
    open(file) do f
        readline(f) # skip header
        for line in eachline(f)
            fields = split(line, '\t')
            template = FeatureTemplate(fields[1], parse(Bool, fields[2]), parse(Float32, fields[3]))
            templates[template.path] = template
        end
    end
    return templates
end

struct FeatureStack
    path::String # == template.path
    stack::CircularVector
    template::FeatureTemplate
end
datasize(f::FeatureStack) = sizeof(FeatureStack) + sizeof(f.path) + length(f.stack) * sizeof(Int8)

DFeatureStack = Dict{String,FeatureStack}
AFeatureStack = Vector{FeatureStack}

function fill_feature_stack(target_length::Int32, annotations::Vector{Annotation},
    feature_templates::Dict{String,FeatureTemplate})::AFeatureStack
    # assumes annotations are ordered by path
    # so that each stacks[idx] array has the same annotation.path
    stacks = AFeatureStack()
    sizehint!(stacks, length(feature_templates))
    indexmap = Dict{String,Int}()
    for annotation in annotations
        template = get(feature_templates, annotation.path, nothing)
        if template === nothing
            @error "Can't find template for $(annotation.path)"
            continue
        end
        index = get(indexmap, annotation.path, nothing)
        if isnothing(index)
            stack = FeatureStack(annotation.path, CircularVector(zeros(Int8, target_length)), template)
            push!(stacks, stack)
            indexmap[annotation.path] = length(stacks)
        else
            stack = stacks[index]
        end
        stack_stack = stack.stack

        @inbounds for i = annotation.start:annotation.start+annotation.length-one(Int32)
            stack_stack[i] += one(Int8)
        end
    end
    @debug "found ($(length(stacks))) $(human(datasize(stacks))) FeatureStacks from $(length(annotations)) annotations"
    return stacks
end

function align_template(feature_stack::FeatureStack)::Tuple{Vector{Tuple{Int32,Int}},Int32}
    stack = feature_stack.stack
    glen = length(stack)
    median_length = floor(Int32, feature_stack.template.median_length)

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

function refine_boundaries_by_score!(feat1::Feature, feat2::Feature, genome_length::Int32, stacks::DFeatureStack)
    # feat1 should be before feat2
    start = feat1.start + feat1.length - one(Int32)
    range_to_test = min(start, feat2.start):max(start, feat2.start)

    feat1_stack = stacks[annotation_path(feat1)].stack
    feat2_stack = stacks[annotation_path(feat2)].stack

    score = sum(feat2_stack[range_to_test]) - sum(feat1_stack[range_to_test])
    maxscore = score
    fulcrum = first(range_to_test) - one(Int32)
    @inbounds for i in range_to_test
        score += 2 * feat1_stack[i]
        score -= 2 * feat2_stack[i]
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

function refine_gene_models!(gene_models::Vector{Vector{SFF_Feature}}, target_seq::CircularSequence,
    feature_stacks::DFeatureStack, orfs::Vector{Feature})

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
                    refine_boundaries_by_score!(previous_feature, feature, genome_length, feature_stacks)
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

mutable struct SFF_Model
    gene::String
    gene_prob::Float32 # mean of probs of component features
    strand::Char
    gene_count::Int8
    exon_count::Int8
    features::Vector{SFF_Feature}
    warnings::Vector{String}
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

gene_length(model::Vector{Feature}) = begin
    maximum([f.start + f.length for f in model]) - minimum([f.start for f in model])
end

gene_length(model::SFF_Model) = begin
    maximum([m.feature.start + m.feature.length for m in model]) - minimum([m.feature.start for m in model])
end

const LACKS_ESSENTIAL_FEATURE = "lacks an essential part of the gene"
const LACKS_START_CODON = "lacks a start codon"
const PREMATURE_STOP_CODON = "has a premature stop codon"
const CDS_NOT_DIVISIBLE_BY_3 = "CDS is not divisible by 3"
const LOW_ANNOTATION_DEPTH = "has low annotation depth, probably spurious"
const INFERIOR_COPY = "probable pseudogene as better copy exists in the genome"
const PARTIAL_IR = "probable pseudogene as overlaps end of inverted repeat"
const OVERLAPPING_FEATURE = "better-scoring feature overlaps with this one"

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
function feature_xgb(ftype::String, template::FeatureTemplate, featurelength::Int32, fdepth::Float32, gmatch::Float32, codingprob::Float32)::Float32
    featurelength ≤ 0 && return Float32(0.0)
    fdepth ≤ 0 && return Float32(0.0)
    rtlength = template.median_length / MAXFEATURELENGTH
    rflength = featurelength / MAXFEATURELENGTH
    frtlength = featurelength / template.median_length
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

function reverse_complement(r::UnitRange, glength::Int32)::UnitRange
    return range(mod1(glength - r.stop + 1, glength), length=length(r))
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
                elseif featuretype(model1.features) == "CDS"  # overlap IR, keep longest model
                    if length(model1_boundaries) > length(model2_boundaries)
                        push!(model2.warnings, PARTIAL_IR)
                    elseif length(model2_boundaries) > length(model1_boundaries)
                        push!(model1.warnings, PARTIAL_IR)
                    end
                else                                           # overlap IR, keep best model
                    if model1.gene_prob > model2.gene_prob
                        push!(model2.warnings, PARTIAL_IR)
                    elseif model2.gene_prob > model1.gene_prob
                        push!(model1.warnings, PARTIAL_IR)
                    end
                end
            end
        end
    end
end

MaybeIR = Union{AlignedBlock,Nothing}

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
    finish %= genome_length

    if strand == '-'
        start = genome_length - finish + 1
        finish = start + length - 1
        start %= genome_length
        finish %= genome_length
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

include("rnas.jl")

