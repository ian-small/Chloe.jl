import Printf: @sprintf
import StatsBase
import IntervalTrees: IntervalBTree, AbstractInterval, Interval
using BioSequences
using JLD2

mutable struct Feature <: AbstractInterval{Int32}
    gene::String
    type::String
    order::Char
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
        this.order = length(tags) == 4 ? tags[4][1] : tags[3][1]
        this.start = s
        this.length = l
        this.phase = p
        this
    end
end

# number of features we need to scan
RefTree = IntervalBTree{Int32,Feature,64}

Base.first(f::Feature) = f.start
function Base.last(f::Feature)
    f.start + f.length - one
end
Base.intersect(i::RefTree, lo::Integer, hi::Integer) = intersect(i, Interval{Int32}(lo, hi))

datasize(f::Feature) = sizeof(Feature) + sizeof(f.annotationPath()) + sum(sizeof(p) for p in f._path_components)
datasize(i::RefTree) = sum(datasize(f) for f in i)

const annotationPath(f::Feature) = begin
    join([f.gene, f.type, string(f.order)], "/")
end

AFeature = Vector{Feature}
AAFeature = Vector{AFeature}
# entire set of Features for one strand of one genome
struct FeatureArray
    genome_id::String
    genome_length::Int32
    strand::Char
    # features::AFeature
    interval_tree::RefTree
end
datasize(f::FeatureArray) = begin
    sizeof(FeatureArray) + sizeof(f.genome_id)  + datasize(f.interval_tree)
end

#extended Feature struct to add prediction info
mutable struct SFF_Feature
    feature::Feature
    relative_length::Float32
    stackdepth::Float32
    gmatch::Float32
    feature_prob::Float32
    coding_prob::Float32
end

function readFeatures(file::String)::FwdRev{FeatureArray}
    open(file) do f
        header = split(readline(f), '\t')
        genome_id = header[1]
        genome_length = parse(Int32, header[2])
        r_features = AFeature()
        f_features = AFeature()
        while !eof(f)
            fields = split(readline(f), '\t')
            startswith(fields[1], "unassigned") && continue #don't use unassigned annotations
            startswith(fields[1], "predicted") && continue #don't use annotations considered as predictions
            length(fields) ≥ 9 && occursin("pseudo", fields[9]) && continue #don't use annotations considered as pseudogenes
            feature = Feature(fields[1], parse(Int, fields[3]), parse(Int, fields[4]), parse(Int, fields[5]))
            if fields[2][1] == '+'
                push!(f_features, feature)
            else
                push!(r_features, feature)
            end
        end
        sort!(f_features)
        sort!(r_features)
        f_strand_features = FeatureArray(genome_id, genome_length, '+', RefTree(f_features))
        r_strand_features = FeatureArray(genome_id, genome_length, '-', RefTree(r_features))
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
datasize(a::Annotation) = sizeof(Annotation)  + sizeof(a.path)# genome_id is shared + sizeof(a.genome_id)
datasize(v::Vector{Annotation}) = sum(datasize(a) for a in v)

function addAnnotation(genome_id::String, feature::Feature, block::AlignedBlock)::Union{Nothing,Annotation}
    #println(block)
    #coordinates in source genome
    src_start = max(feature.start, block.src_index)
    src_end::Int32 = min(feature.start + feature.length, block.src_index + block.blocklength) - 1
    src_length::Int32 = src_end - src_start + 1
    src_length <= 0 && return nothing
    offset5 = src_start - feature.start
    offset3 = src_end - (feature.start + feature.length -1)
    if (feature.type == "CDS")
        phase = phaseCounter(feature.phase, offset5)
    else
        phase = 0
    end
    #coordinates in target genome
    tgt_start = src_start > feature.start ? block.tgt_index : block.tgt_index + (feature.start - block.src_index)

    annotation_path = annotationPath(feature)
    Annotation(genome_id, annotation_path, 
                            tgt_start,
                            src_length, offset5,
                            offset3,
                            phase)
end

function findOverlaps(ref_features::FeatureArray, aligned_blocks::AlignedBlocks)::Vector{Annotation}
    annotations = Vector{Annotation}()
    for block in aligned_blocks
        for feature in intersect(ref_features.interval_tree, block.src_index, block.src_index + block.blocklength - one)
            # @assert rangesOverlap(feature.start, feature.length, block[1], block[3])
            anno = addAnnotation(ref_features.genome_id, feature, block)
            if anno !== nothing
                push!(annotations, anno)
            end
        end
    end
    annotations
end

struct FeatureTemplate
    path::String  # similar to .sff path
    median_length::Float32 #median length of feature
    #GLM coefficients
    intercept::Float64
    length_coef::Float64
    depth_coef::Float64
    match_coef::Float64
    depth_length_coef::Float64
    depth_match_coef::Float64
end

datasize(s::FeatureTemplate) = sizeof(FeatureTemplate) + sizeof(s.path)


function readTemplates(file::String)::Dict{String,FeatureTemplate}

    if !isfile(file)
        error("\"$(file)\" is not a file")
    end
    if filesize(file) === 0
        error("no data in \"$(file)!\"")
    end

    gene_exons = String[]
    # templates = FeatureTemplate[]
    templates = Dict{String,FeatureTemplate}()

    open(file) do f
        header = readline(f)
        for line in eachline(f)
            fields = split(line, '\t')
            path_components = split(fields[1], '/')
            if path_components[2] ≠ "intron"
                push!(gene_exons, path_components[1])
            end
            template = FeatureTemplate(fields[1], parse(Float32, fields[2]), parse(Float64, fields[3]), parse(Float64, fields[4]),
                 parse(Float64, fields[5]), parse(Float64, fields[6]), parse(Float64, fields[7]), parse(Float64, fields[8]))
            if haskey(templates, template.path)
                @error "duplicate path: $(template.path) in \"$file\""
            end
            templates[template.path] = template
            # push!(templates, template)
        end
    end
    # sort!(templates, by=x -> x.path)

    return templates
end

struct FeatureStack
    path::String # == template.path
    stack::CircularVector
    template::FeatureTemplate
end
datasize(f::FeatureStack) = sizeof(FeatureStack) + sizeof(f.path) + length(f.stack) * sizeof(Int32)

DFeatureStack = Dict{String,FeatureStack}
AFeatureStack = Vector{FeatureStack}
# AFeatureStack = Dict{String,FeatureStack}
ShadowStack = CircularVector

function fillFeatureStack(target_length::Int32, annotations::Vector{Annotation},
    feature_templates::Dict{String,FeatureTemplate})::Tuple{AFeatureStack,ShadowStack}

    # Implementation Note: Feature stacks are rather large (we can get 80MB)
    # ... but! the memory requirement is upperbounded by the number of features
    # in the optimized_templates.v2.tsv

    # assumes annotations is ordered by path
    # so that each stacks[idx] array has the same annotation.path
    stacks = AFeatureStack()
    sizehint!(stacks, length(feature_templates))
    indexmap = Dict{String,Int}()
    shadowstack::ShadowStack = ShadowStack(zeros(Int32, target_length)) # will be negative image of all stacks combined,
    for annotation in annotations
        template = get(feature_templates, annotation.path, nothing)
        # template_index = findfirst(x -> x.path == annotation.path, templates)
        if template === nothing
            #@error "Can't find template for $(annotation.path)"
            continue
        end
        index = get(indexmap, annotation.path, nothing)
        if index === nothing
            stack = FeatureStack(annotation.path, CircularVector(zeros(Int32, target_length)), template)
            push!(stacks, stack)
            indexmap[annotation.path] = length(stacks)
        else
            stack = stacks[index]
        end
        stack_stack = stack.stack

        @inbounds for i = annotation.start:annotation.start + annotation.length - one
            stack_stack[i] += one
            shadowstack[i] -= one
        end
    end
    @debug "found ($(length(stacks))) $(human(datasize(stacks) + datasize(shadowstack))) FeatureStacks from $(length(annotations)) annotations"
    return stacks, shadowstack
end

function expandBoundaryInChunks(feature_stack::FeatureStack, shadowstack::ShadowStack, origin, direction, max)::Int
    glen = length(shadowstack)
    index = origin

    # expand in chunks
    for chunksize in [100,80,60,40,30,20,10,5,1]
        for i = direction:direction * chunksize:direction * max
            sumscore = 0
            for j = 1:chunksize
                sumscore += feature_stack.stack[index + direction * j] + shadowstack[index + direction * j]
            end
            sumscore < 0 && break
            index += direction * chunksize
        end
    end
    return genome_wrap(glen, index)
end

function expandBoundary(feature_stack::FeatureStack, shadowstack::ShadowStack, origin, direction, max)
    glen = length(shadowstack)
    stack_stack = feature_stack.stack
    index = origin
    @inbounds for i = direction:direction:direction * max
        sumscore = stack_stack[index + direction] + shadowstack[index + direction]
        sumscore < 0 && break
        index += direction
    end
    return index # not wrapped
end

function getDepthAndCoverage(feature_stack::FeatureStack, left::Int32, len::Int32, numrefs::Int)::Tuple{Float64,Float64}
    coverage = 0
    max_count = 0
    sum_count = 0
    stack_stack = feature_stack.stack
    stack_len = length(stack_stack)
    @inbounds for nt::Int32 = left:left + len - one
        count = stack_stack[nt]
        if count > 0
            coverage += 1
            sum_count += count
        end
        if count > max_count
            max_count = count
        end
    end
    if max_count == 0
        depth = 0
    else
        depth = sum_count / (numrefs * len)
    end
    return depth, coverage / len
end

function alignTemplateToStack(feature_stack::FeatureStack)::Tuple{Vector{Tuple{Int32,Int}}, Int32}
    stack = feature_stack.stack
    glen = length(stack)
    median_length = floor(Int32, feature_stack.template.median_length)

    hits = Vector{Tuple{Int32,Int}}()
    score = 0
    @inbounds for nt::Int32 = one:median_length
        score += stack[nt]
    end

    if score > 0; push!(hits, (one, score)); end
    
    @inbounds for nt::Int32 = 2:glen
        mt::Int32 = nt + median_length - one
        score -= stack[nt - one] # remove tail
        score += stack[mt] # add head
        if score > 0
            if isempty(hits)
                push!(hits, (nt, score))
            else
                previous = last(hits)
                if nt ≥ previous[1] + median_length
                    push!(hits, (nt, score))
                elseif score > previous[2]  #if hits overlap, keep the best one
                    pop!(hits)
                    push!(hits, (nt, score))
                end
            end   
        end
    end
    isempty(hits) && return hits, median_length
    #sort by descending score
    sort!(hits, by = x -> x[2], rev = true)
    maxscore = hits[1][2]
    #retain all hits scoring over 90% of the the top score
    filter!(x -> x[2] ≥ maxscore*0.9, hits)
    return hits, median_length
end

#= function getFeaturePhaseFromAnnotationOffsets(feat::Feature, annotations::Vector{Annotation})::Int8
    phases = Int8[]
    for annotation in annotations[findall(x -> x.path == annotationPath(feat), annotations)]
        if rangesOverlap(feat.start, feat.length, annotation.start, annotation.length)
            # estimating feature phase from annotation phase
            phase = phaseCounter(annotation.phase, feat.start - annotation.start)
            push!(phases, phase)
        end
    end
    if length(phases) == 0
        return Int8(0)
    end
    println(feat, '\t', phases)
    return StatsBase.mode(phases) # return most common phase
end =#

function weightedMode(values::Vector{Int32}, weights::Vector{Float32})::Int32
    w, v = findmax(StatsBase.addcounts!(Dict{Int32,Float32}(), values, weights))
    return v
end

# uses weighted mode, weighting by alignment length and distance from boundary
function refineMatchBoundariesByOffsets!(feat::Feature, annotations::Vector{Annotation}, 
            target_length::Integer, coverages::Dict{String,Float32})
    # grab all the matching features and sort by start in genome of origin to ensure that when iterated in order, most 5' match is found first
    matching_annotations = sort(annotations[findall(x -> x.path == annotationPath(feat), annotations)], by = a -> a.start)
    isempty(matching_annotations) && return #  feat, [], []
    overlapping_annotations = Annotation[]
    minstart = target_length
    maxend = 1

    ## Fix 5' boundary and feature phase
    for annotation in matching_annotations
        annotation.offset5 >= feat.length && continue
        if rangesOverlap(feat.start, feat.length, annotation.start, annotation.length)
            #ignore if we already have an annotation from this genome
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
        end5ws[i] =  weight * weight
        phases[i] = Int32(phaseCounter(annotation.phase, -annotation.offset5))
        #println(feat, '\t', annotation)
    end
    left = feat.start
    if !isempty(end5s)
        left = weightedMode(end5s, end5ws)
    end
    feat.start = genome_wrap(target_length, left)
    feat.phase = isempty(phases) ? 0 : weightedMode(phases, end5ws)
    #println(annotationPath(feat) , '\t', phases, '\t', feat.phase)

    ## Fix 3' boundary
    empty!(overlapping_annotations)
    for annotation in Iterators.reverse(matching_annotations) #iterate in reverse to ensure 3' annotations are reached first
        annotation.offset3 >= feat.length && continue
        if rangesOverlap(feat.start, feat.length, annotation.start, annotation.length)
            #ignore if we already have an annotation from this genome
            if isnothing(findfirst(a -> a.genome_id == annotation.genome_id, overlapping_annotations))
                maxend = max(maxend, annotation.start + annotation.length -1)
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
        right = weightedMode(end3s, end3ws)
    end
    feat.length = right - left + 1
end

function groupFeaturesIntoGeneModels(sff_features::Vector{SFF_Feature})::Vector{Vector{SFF_Feature}}
    gene_models = [SFF_Feature[]]
    current_model = SFF_Feature[]
    left_border::Int32 = typemax(Int32)
    right_border::Int32 = 0
    for sff_feature in sff_features
        if isempty(current_model)
            push!(current_model, sff_feature)
            left_border = sff_feature.feature.start
            right_border = sff_feature.feature.start + sff_feature.feature.length - 1
        #add feature to model if the gene names match and they are less than 3kb apart; 3kb is an arbitrary limit set to include the longest intron (trnK, ~2.5kb) and may need to be reduced
        #a shorter distance limit would require the features to be pre-sorted by position, not name, so that we know they are being added in order
        elseif current_model[1].feature.gene == sff_feature.feature.gene && min(left_border - (sff_feature.feature.start + sff_feature.feature.length - 1), sff_feature.feature.start - right_border) < 3000
            push!(current_model, sff_feature)
            left_border = min(left_border, sff_feature.feature.start)
            right_border = max(right_border, sff_feature.feature.start + sff_feature.feature.length - 1)
        else
            sort!(current_model, by = x -> x.feature.start)
            push!(gene_models, current_model)
            current_model = SFF_Feature[]
            push!(current_model, sff_feature)
            left_border = sff_feature.feature.start
            right_border = sff_feature.feature.start + sff_feature.feature.length - 1
        end
    end
    if length(current_model) > 0
        sort!(current_model, by = x -> x.feature.start)
        push!(gene_models, current_model)
    end
    return gene_models
end

function translateModel(target_seq::CircularSequence, model::Vector{SFF_Feature})::LongAminoAcidSeq
    cds = dna""
    empty!(cds)
    for sfeat in model
        feat = sfeat.feature
        feat.type ≠ "CDS" && continue
        append!(cds, target_seq[feat.start:feat.start + feat.length - 1])
    end
    cds = cds[first(model).feature.phase + 1:end]
    if mod(length(cds), 3) ≠ 0
        println(model)
    end
    return translate(cds)
end

# define parameter weights
# const start_score = Dict("ATG"=>1.0,"ACG"=>0.1,"GTG"

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

# function findStartCodon2!(cds::Feature, genome_length::Integer, genomeloop::CircularSequence, predicted_starts,
#                           predicted_starts_weights, feature_stack, shadowstack)
#     # define range in which to search from upstream stop to end of feat
#     start5 = genome_wrap(genome_length, cds.start + cds.phase)
#     codon = SubString(genomeloop, start5, start5 + 2)
#     while !isStopCodon(codon, false)
#         start5 = genome_wrap(genome_length, start5 - 3)
#         codon = SubString(genomeloop, start5, start5 + 2)
#     end
#     search_range = start5:cds.start + cds.length - 1
#     window_size = length(search_range)

#     codons = zeros(window_size)
#     phase = zeros(window_size)
#     cumulative_stack_coverage = zeros(window_size)
#     cumulative_shadow_coverage = zeros(window_size)

#     stack_coverage = 0
#     shadow_coverage = 0
#     for (i, nt) in enumerate(search_range)
#         gw = genome_wrap(genome_length, nt)
#         codon = SubString(genomeloop, nt, nt + 2)
#         codons[i] = startScore(cds, codon)
#         if i % 3 == 1
#             phase[i] = 1.0
#         end
#         if feature_stack[gw] > 0
#             stack_coverage += 1
#         elseif shadowstack[gw] < 0
#             shadow_coverage -= 1
#         end
#         if stack_coverage > 3
#             cumulative_stack_coverage[i] = 3 - stack_coverage
#         end
#         cumulative_shadow_coverage[i] = shadow_coverage
#     end

#     predicted_starts = zeros(window_size)
#     for (s, w) in zip(predicted_starts, predicted_starts_weights)
#         for (i, nt) in enumerate(search_range)
#             distance = 1 + (s - nt)^2
#             predicted_starts[i] = w / distance
#         end
#     end

#     # combine vectors using parameter weights
#     result = codons .* phase .* cumulative_stack_coverage .+ cumulative_shadow_coverage .+ predicted_starts
#     # set feat.start to highest scoring position
#     maxstartscore, maxstartpos = findmax(result)
#     maxstartpos += start5 - 1
#     cds.length += cds.start - maxstartpos
#     cds.start = genome_wrap(genome_length, maxstartpos)
#     # set phase to zero
#     cds.phase = 0
#     return cds
# end

function findStartCodon!(cds::Feature, genome_length::Int32, target_seq::CircularSequence)
    # assumes phase has been correctly set
    # search for start codon 5'-3' beginning at cds.start, save result; abort if stop encountered
    start3::Int32 = cds.start + cds.phase
    codon = getcodon(target_seq, start3)
    if isStartCodon(codon, true, true)  # allow ACG and GTG codons if this is predicted start
        cds.start = start3
        cds.length -= cds.phase
        cds.phase = 0;
        return cds
    end

    allowACG = false
    allowGTG = false
    if cds.gene == "rps19"
        allowGTG = true
    end
    
    start3 = 0
    for i::Int32 in cds.start + cds.phase:3:cds.start + cds.phase + genome_length - 3
        codon = getcodon(target_seq, i)
        if isStartCodon(codon, allowACG, allowGTG)
            start3 = i >= cds.start + cds.phase ? i : i + genome_length
            break
        end
        if isStopCodon(codon, false)
            break
        end
    end
    # search for start codon 3'-5' beginning at cds.start, save result; abort if stop encountered
    start5::Int32 = 0
    for i::Int32 in cds.start + cds.phase:-3:cds.start + cds.phase - genome_length + 3
        codon = getcodon(target_seq, i)
        if isStartCodon(codon, allowACG, allowGTG)
            start5 = i <= cds.start + cds.phase ? i : i - genome_length
            break
        end
        if isStopCodon(codon, false)
            break
        end
    end
    # return cds with start set to nearer of the two choices
    if start3 == 0 && start5 == 0
        #Couldn't find start
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

function setLongestORF!(sfeat::SFF_Feature, orfs::Vector{Feature})
    #find in-frame orf with longest overlap with feature
    #orfs are sorted by descending length
    #println(sfeat)
    f = sfeat.feature
    original_start = f.start
    original_length = f.length
    f_frame = mod1(f.start + f.phase, 3)
    max_overlap = 0
    overlap_orf = nothing
    for orf in orfs
        orf.length < max_overlap && break
        orf_frame = mod1(orf.start, 3)
        if f_frame == orf_frame && rangesOverlap(orf.start, orf.length, f.start, f.length)
            overlap = length(intersect(range(orf.start, length = orf.length), range(f.start, length = f.length)))
            if overlap > max_overlap
                max_overlap = overlap
                overlap_orf = orf
            end
        end
    end
    isnothing(overlap_orf) && return
    if overlap_orf.start > f.start # must be an internal stop
        f.phase = 0
        f.start = overlap_orf.start
    end
    f.length = overlap_orf.start - f.start + overlap_orf.length
end

function refineBoundariesbyScore!(feat1::Feature, feat2::Feature, genome_length::Int32, stacks::DFeatureStack)
    # feat1 should be before feat2
    start = feat1.start + feat1.length - one
    range_to_test = min(start, feat2.start):max(start, feat2.start)
    
    feat1_stack = stacks[annotationPath(feat1)].stack
    feat2_stack = stacks[annotationPath(feat2)].stack

    # score = sum(@view feat2_stack[range_to_test]) - sum(@view feat1_stack[range_to_test])
    score = sum(feat2_stack[range_to_test]) - sum(feat1_stack[range_to_test])
    maxscore = score
    fulcrum = first(range_to_test) - one
    @inbounds for i in range_to_test
        score += 2 * feat1_stack[i]
        score -= 2 * feat2_stack[i]
        if score > maxscore
            maxscore = score
            fulcrum = i
        end
    end
    # fulcrum should point to end of feat1
    feat1.length = fulcrum - feat1.start + 1
    # fulcrum+1 should point to start of feat2
    feat2.length += feat2.start - (fulcrum + 1)
    feat2.start = genome_wrap(genome_length, fulcrum + 1)
end

function refineGeneModels!(gene_models::Vector{Vector{SFF_Feature}}, target_seq::CircularSequence,
                          annotations::Vector{Annotation},
                          feature_stacks::DFeatureStack,
                          orfs::Vector{Feature})

    ## add code to make use of the feature_prob field in the SFF_Features

    genome_length = Int32(length(target_seq))
    last_cds_examined = nothing
    for model in gene_models
        isempty(model) && continue
        # sort features in model by mid-point to avoid cases where long intron overlaps short exon
        sort!(model, by=f -> f.feature.start + f.feature.length / 2)
        last_exon = last(model).feature
        # if CDS, find phase, stop codon and set feature.length
        if last_exon.type == "CDS"
            #last_exon.phase = getFeaturePhaseFromAnnotationOffsets(last_exon, annotations)
            if last_exon.gene ≠ "rps12A"
                setLongestORF!(last(model), orfs)
            end
            last_cds_examined = last_exon
        end
        features_to_remove = SFF_Feature[]
        for i in length(model) - 1:-1:1
            feature = model[i].feature
            next_feature = model[i + 1].feature
            # check adjacent to last exon, if not...
            gap = next_feature.start - (feature.start + feature.length)
            if gap ≠ 0 && gap < 100
                refineBoundariesbyScore!(feature, next_feature, genome_length, feature_stacks)
                gap = next_feature.start - (feature.start + feature.length)
            end
            feature.length ≤ 0 && push!(features_to_remove, model[i])
            next_feature.length ≤ 0 && push!(features_to_remove, model[i+1])
            # if CDS, check phase is compatible
            if feature.type == "CDS"
                #feature.phase = getFeaturePhaseFromAnnotationOffsets(feature, annotations)
                if last_cds_examined !== nothing && phaseCounter(feature.phase, feature.length % Int32(3)) != last_cds_examined.phase
                    # refine boundaries (how?)
                end
                last_cds_examined = feature
            end
        end
        setdiff!(model, features_to_remove)
        # if CDS, find start codon and set feature.start
        first_exon = first(model).feature
        if first_exon.type == "CDS"
            #first_exon.phase = getFeaturePhaseFromAnnotationOffsets(first_exon, annotations)
            if last_exon.gene ≠ "rps12B"
                findStartCodon!(first_exon, genome_length, target_seq)
                # first_exon = findStartCodon2!(first_exon,genome_length,targetloop)
            end
        end
    end
end

mutable struct SFF_Model
    gene::String
    gene_prob::Float32 #mean of probs of comonent features
    strand::Char
    gene_count::Int8
    exon_count::Int8
    features::Vector{SFF_Feature}
    hasStart::Bool
    hasPrematureStop::Bool
end

function get_gene_boundaries(model::SFF_Model)::UnitRange{Int32}
    start = first(model.features).feature.start
    length = last(model.features).feature.start + last(model.features).feature.length - start
    return range(start, length = length)
end

function mean_stackdepth(model::SFF_Model)::Float64
    sum = 0.0
    for sff in model.features
        sum += sff.stackdepth
    end
    return sum/length(model.features)
end

gene_length(model::Vector{Feature}) = begin
    maximum([f.start + f.length for f in model]) - minimum([f.start for f in model])
end

gene_length(model::Vector{SFF_Model}) = begin
    maximum([m.feature.start + m.feature.length for m in model]) - minimum([m.feature.start for m in model])
end

function toSFF(model::Vector{SFF_Feature}, strand::Char, target_seq::CircularSequence, sensitivity::Float16)::Union{Nothing, SFF_Model}
    first_model = first(model)
    gene = first_model.feature.gene
    exon_count = 0
    cds = false
    model_prob = 0.0
    for (n, sff_feature) in enumerate(model)
        #keep model if mean feature probability exceeds the sensitivity threshold
        model_prob = n == 1 ? 1 - sff_feature.feature_prob : model_prob * (1 - sff_feature.feature_prob)
        if sff_feature.feature.type ≠ "intron"
            exon_count += 1
        end
    end
    model_prob = 1 - model_prob
    exceeds_sensitivity = false
    if model_prob ≥ sensitivity || isnan(model_prob)
        exceeds_sensitivity = true
    end
    hasStart = true
    hasPrematureStop = false
    if exceeds_sensitivity && cds
        #protein = translateModel(target_seq, model)
        #println(gene, '\t', protein)
        if gene ≠ "rps12B" && !isStartCodon(getcodon(target_seq, first_model.feature.start), true, true)
            hasStart = false
        end
        #stop_position = findfirst(AA_Term, protein)
        #println(stop_position)
        #if !isnothing(stop_position) && stop_position < length(protein)
        #    hasPrematureStop = true
        #end
    end
    return exceeds_sensitivity ? SFF_Model(gene, model_prob, strand, 1, exon_count, model, hasStart, hasPrematureStop) : nothing
end

function writeModelToSFF(outfile::IO, model::SFF_Model)

    isnothing(model) && return
    model_id = model.gene * "/" * string(model.gene_count)
    !model.hasStart && !startswith(model.gene, "IR") && @warn("$(model_id) has no start codon")
    model.hasPrematureStop && @warn("$(model_id) has (a) premature stop codon(s)")
    
    for sff in model.features
        f = sff.feature
        write(outfile, model_id)
        write(outfile, "/")
        write(outfile, f.type)
        write(outfile, "/")
        write(outfile, string(f.order))
        write(outfile, "\t")
        write(outfile, join([model.strand,string(f.start),string(f.length),string(f.phase)], "\t"))
        write(outfile, "\t")
        write(outfile, join([@sprintf("%.3f",sff.relative_length),@sprintf("%.3f",sff.stackdepth),@sprintf("%.3f",sff.gmatch),
            @sprintf("%.3f",sff.feature_prob), @sprintf("%.3f",sff.coding_prob)], "\t"))
        write(outfile, "\n")
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

function feature_glm(template::FeatureTemplate, flength::Float32, fdepth::Float32, gmatch::Float32)::Float32
    pred = template.intercept + template.length_coef * flength + template.depth_coef * fdepth + template.match_coef * gmatch + fdepth * flength * template.depth_length_coef + fdepth * gmatch * template.depth_match_coef
    odds = exp(pred)
    return Float32(odds / (1.0 + odds))
end

# only some gene_models are allowed to overlap
function allowed_model_overlap(m1, m2)::Bool
    if m1.gene == "trnK-UUU" && m2.gene == "matK"; return true; end
    if m1.gene == "matK" && m2.gene == "trnK-UUU"; return true; end
    if m1.gene == "ndhC" && m2.gene == "ndhK"; return true; end
    if m1.gene == "ndhK" && m2.gene == "ndhC"; return true; end
    if m1.gene == "psbC" && m2.gene == "psbD"; return true; end
    if m1.gene == "psbD" && m2.gene == "psbC"; return true; end
    if m1.gene == "atpB" && m2.gene == "atpE"; return true; end
    if m1.gene == "atpE" && m2.gene == "atpB"; return true; end
    if m1.gene == "rps3" && m2.gene == "rpl22"; return true; end
    if m1.gene == "rpl22" && m2.gene == "rps3"; return true; end
    return false
end

function filter_gene_models!(fwd_models::Vector{SFF_Model}, rev_models::Vector{SFF_Model}, numrefs::Int)
    #filter out models of genes where better model of same gene exists
    #filter out models of genes where better model of overlapping gene exists
    #filter on stackdepth(sanity check in case of GLM failure)

    fwd_models_to_delete = SFF_Model[]
    for (i, model1) in enumerate(fwd_models)
        if mean_stackdepth(model1) < 1/(4 * numrefs); push!(fwd_models_to_delete, model1); continue; end
        for j in i + 1:length(fwd_models)
            model2 = fwd_models[j]
            bmodel1 = get_gene_boundaries(model1)
            bmodel2 = get_gene_boundaries(model2)
            if model1.gene == model2.gene || (length(intersect(bmodel1, bmodel2))> 0 && !allowed_model_overlap(model1, model2))
                if model1.gene_prob > 1.05 * model2.gene_prob ##small tolerance for unequal probs
                    push!(fwd_models_to_delete, model2)
                elseif model2.gene_prob > 1.05 * model1.gene_prob
                    push!(fwd_models_to_delete, model1)
                end
            end
        end
    end
    #println(fwd_models_to_delete)
    setdiff!(fwd_models, fwd_models_to_delete)
    
    rev_models_to_delete = SFF_Model[]
    for (i, model1) in enumerate(rev_models)
        if mean_stackdepth(model1) < 1/(4 * numrefs); push!(rev_models_to_delete, model1); continue; end
        for j in i + 1:length(fwd_models)
            model2 = fwd_models[j]
            bmodel1 = get_gene_boundaries(model1)
            bmodel2 = get_gene_boundaries(model2)
            if model1.gene == model2.gene || (length(intersect(bmodel1, bmodel2)) > 0 && bmodel1[1] % 3 == bmodel2[1] % 3)
                if model1.gene_prob > 1.05 * model2.gene_prob ##small tolerance for unequal probs
                    push!(rev_models_to_delete, model2)
                elseif model2.gene_prob > 1.05 * model1.gene_prob
                    push!(rev_models_to_delete, model1)
                end
            end
        end
    end
    setdiff!(rev_models, rev_models_to_delete)

    empty!(fwd_models_to_delete)
    empty!(rev_models_to_delete)
    for model1 in fwd_models, model2 in rev_models
        if model1.gene == model2.gene
            if model1.gene_prob > 1.05 * model2.gene_prob ##small tolerance for unequal probs
                push!(rev_models_to_delete, model2)
            elseif model2.gene_prob > 1.05 * model1.gene_prob
                push!(fwd_models_to_delete, model1)
            end
        end
    end
    println(fwd_models_to_delete)
    setdiff!(fwd_models, fwd_models_to_delete)
    setdiff!(rev_models, rev_models_to_delete)
end

MaybeIR = Union{AlignedBlock,Nothing}

function getModelID!(model_ids::Dict{String,Int32}, model::SFF_Model)
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

function writeSFF(outfile::Union{String,IO},
                id::String, # NCBI id
                genome_length::Int32,
                mean_coverage::Float32,
                models::FwdRev{Vector{SFF_Model}})

    function out(outfile::IO)
        model_ids = Dict{String,Int32}()
        write(outfile, id, "\t", string(genome_length), "\t", @sprintf("%.3f",mean_coverage), "\n")
        for model in models.forward
            isnothing(model) && continue
            isempty(model.features) && continue
            model_id = getModelID!(model_ids, model)
            writeModelToSFF(outfile, model)
        end
        for model in models.reverse
            isnothing(model) && continue
            isempty(model.features) && continue
            model_id = getModelID!(model_ids, model)
            writeModelToSFF(outfile, model)
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
    
    if strand == '-'
        start = genome_length - finish + 1
        if start < 1
            start += genome_length
        end
        finish = start + length - 1
    end
    return (start, finish, length)
end

function writeGFF3(outfile::Union{String,IO},
    genome_id::String, # NCBI id
    genome_length::Int32,
    models::FwdRev{Vector{SFF_Model}})

    function out(outfile::IO)
        model_ids = Dict{String,Int32}()
        write(outfile, "##gff-version 3.2.1\n")
        for model in models.forward
            isnothing(model) && continue
            isempty(model.features) && continue
            model.gene_count = getModelID!(model_ids, model)
        end
        for model in models.reverse
            isnothing(model) && continue
            isempty(model.features) && continue
            model.gene_count = getModelID!(model_ids, model)
        end
        allmodels = sort(vcat(models.forward, models.reverse), by = m -> sff2gffcoords(first(m.features).feature, m.strand, genome_length)[1])
        for model in allmodels
            writeModelToGFF3(outfile, model, genome_id, genome_length)
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

function mergeAdjacentFeaturesinModel!(model::SFF_Model)
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

function writeModelToGFF3(outfile, model::SFF_Model, genome_id::String, genome_length::Int32)
    
    function write_line(type, start, finish, pvalue, id, parent, phase="."; key="Parent")
        l = [genome_id, "Chloe", type, start, finish, @sprintf("%.3e",pvalue), model.strand, phase]
        write(outfile, join(l, "\t"))
        write(outfile, "\t", "ID=", id, ";", key, "=", parent, "\n")
    end

    mergeAdjacentFeaturesinModel!(model)

    features = model.features
    id = model.gene
    if model.gene_count > 1
        id = id * "-" * string(model.gene_count)
    end

    start = minimum(f.feature.start for f in features)
    ft = first(features).feature.type
    if ft == "CDS" && !startswith(id, "rps12A")
        last(features).feature.length += 3 #add stop codon
    end
    finish = maximum(f.feature.start + f.feature.length - 1 for f in features)
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
    if ft == "CDS"
        parent =  id * ".mRNA"
        write_line("mRNA", start, finish, 1.0 - model.gene_prob, parent, id)
    elseif ft == "rRNA"
        parent = id * ".rRNA"
        write_line("rRNA", start, finish, 1.0 - model.gene_prob, parent, id)
    elseif ft == "tRNA"
        parent =  id * ".tRNA"
        write_line("tRNA", start, finish, 1.0 - model.gene_prob, parent, id)
    end
    #exons
    for sff in model.features
        f = sff.feature
        type = f.type
        if type == "tRNA" || type == "rRNA"
            type = "exon"
        end
        start, finish, length = sff2gffcoords(f, model.strand, genome_length)
        
        phase = type == "CDS" ? string(f.phase) : "."
        write_line(type, start, finish, 1.0 - sff.feature_prob, annotationPath(f), parent, phase)
    end
    write(outfile, "###\n")
end
