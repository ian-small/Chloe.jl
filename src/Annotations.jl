import Printf: @sprintf
import StatsBase
import IntervalTrees: IntervalBTree, AbstractInterval, Interval

mutable struct Feature <: AbstractInterval{Int32}
    path::String
    start::Int32
    length::Int32
    # phase is the number of nucleotides to skip at the start of the sequence 
    # to be in the correct reading frame
    phase::Int8
    _path_components::Vector{String}
    Feature(path, start, length, phase) = new(path, start, length, phase, split(path, '/'))
end
# number of features we need to scan
RefTree = IntervalBTree{Int32,Feature,64}

Base.first(f::Feature) = f.start
Base.last(f::Feature) = f.start + f.length - one(Int32)
Base.intersect(i::RefTree, lo::Integer, hi::Integer) = intersect(i, Interval{Int32}(lo, hi))

datasize(f::Feature) = sizeof(Feature) + sizeof(f.path) + sum(sizeof(p) for p in f._path_components)
datasize(i::RefTree) = sum(datasize(f) for f in i)

const annotationPath(feature::Feature) = begin
    pc = feature._path_components
    join([pc[1],"?",pc[3], pc[4]], "/")
end
function getFeatureName(feature::Feature)
    feature._path_components[1]
end
# for writing .sff files
function getFeatureSuffix(feature::Feature)
    pc = feature._path_components
    join([pc[3], pc[4]], "/")
end
function getFeatureType(feature::Feature)
    feature._path_components[3]
end

isType(feature::Feature, gene::String) = gene == getFeatureType(feature)

isFeatureName(feature::Feature, name::String) = name == getFeatureName(feature)

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

function readFeatures(file::String)::FwdRev{FeatureArray}
    open(file) do f
        header = split(readline(f), '\t')
        genome_id = header[1]
        genome_length = parse(Int32, header[2])
        r_features = AFeature()
        f_features = AFeature()
        while !eof(f)
            fields = split(readline(f), '\t')
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
    path::String # gene/n/?/n
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
    startA = max(feature.start, block[1])
    flength = min(feature.start + feature.length, block[1] + block[3]) - startA
    flength <= 0 && return nothing
    feature_type = getFeatureType(feature)
    offset5 = startA - feature.start
    if (feature_type == "CDS")
        phase = phaseCounter(feature.phase, offset5)
    else
        phase = 0
    end
    annotation_path = annotationPath(feature)
    Annotation(genome_id, annotation_path, 
                            startA - block[1] + block[2], # shift from reference position to target
                            flength, offset5,
                            feature.length - (startA - feature.start) - flength,
                            phase)
end

# checks all blocks for overlap so could be speeded up by using an interval tree
# function addOverlapBlocks(genome_id::String, feature::Feature, blocks::AlignedBlocks)::Vector{Annotation}
#     annotations = Vector{Annotation}()
#     # sizehint!(pushed_features, length(blocks))

#     for block in blocks
#         if rangesOverlap(feature.start, feature.length, block[1], block[3])
#             anno = addAnnotation(genome_id, feature, block)
#             if anno !== nothing
#                 push!(annotations, anno)
#             end
#         end
#     end
#     annotations
# end

# function findOverlaps(ref_featurearray::FeatureArray, aligned_blocks::AlignedBlocks)::Vector{Annotation}
#     annotations = Vector{Annotation}()
#     for feature in ref_featurearray.features
#         new_annotations = addOverlapBlocks(ref_featurearray.genome_id, feature, aligned_blocks)
#         push!(annotations, new_annotations...)
#     end
#     annotations
# end
function findOverlaps(ref_features::FeatureArray, aligned_blocks::AlignedBlocks)::Vector{Annotation}
    annotations = Vector{Annotation}()
    for block in aligned_blocks
        for feature in intersect(ref_features.interval_tree, block[1], block[1] + block[3] - one(Int32))
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
    threshold_counts::Float32
    threshold_coverage::Float32
    median_length::Int32
end

datasize(s::FeatureTemplate) = sizeof(FeatureTemplate) + sizeof(s.path)


function readTemplates(file::String)::Tuple{Dict{String,FeatureTemplate},Dict{String,Int32}}

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
            if path_components[3] ≠ "intron"
                push!(gene_exons, path_components[1])
            end
            template = FeatureTemplate(fields[1], 
                parse(Float32, fields[2]),
                parse(Float32, fields[3]),
                parse(Int32, fields[4]))
            if haskey(templates, template.path)
                @error "duplicate path: $(template.path) in \"$file\""
            end
            templates[template.path] = template
            # push!(templates, template)
        end
    end
    # sort!(templates, by=x -> x.path)
    return templates, StatsBase.addcounts!(Dict{String,Int32}(), gene_exons)
end

struct FeatureStack
    path::String # == template.path
    stack::Vector{Int32}
    template::FeatureTemplate
end
datasize(f::FeatureStack) = sizeof(FeatureStack) + sizeof(f.path) + length(f.stack) * sizeof(Int32)

DFeatureStack = Dict{String,FeatureStack}
AFeatureStack = Vector{FeatureStack}
# AFeatureStack = Dict{String,FeatureStack}
ShadowStack = Vector{Int32}

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
    shadowstack::ShadowStack = fill(Int32(-1), target_length) # will be negative image of all stacks combined,
    # initialised to small negative number; acts as prior expectation for feature-finding
    for annotation in annotations
        template = get(feature_templates, annotation.path, nothing)
        # template_index = findfirst(x -> x.path == annotation.path, templates)
        if template === nothing
            @error "Can't find template for $(annotation.path)"
            continue
        end
        index = get(indexmap, annotation.path, nothing)
        if index === nothing
            stack = FeatureStack(annotation.path, zeros(Int32, target_length), template)
            push!(stacks, stack)
            indexmap[annotation.path] = length(stacks)
        else
            stack = stacks[index]
        end
        stack_stack = stack.stack
        @inbounds for i = annotation.start:annotation.start + annotation.length - one(Int32)
            gw = genome_wrap(target_length, i)
            stack_stack[gw] += 3 # +1 to counteract shadowstack initialisation, +1 to counteract addition to shadowstack
            shadowstack[gw] -= 1
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
                sumscore += feature_stack.stack[genome_wrap(glen, index + direction * j)] + shadowstack[genome_wrap(glen, index + direction * j)]
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
        gw = genome_wrap(glen, index + direction)
        sumscore = stack_stack[gw] + shadowstack[gw]
        sumscore < 0 && break
        index += direction
    end
    return index # not wrapped
end

function getDepthAndCoverage(feature_stack::FeatureStack, left::Int32, len::Int32)::Tuple{Float64,Float64}
    coverage = 0
    max_count = 0
    sum_count = 0
    stack_stack = feature_stack.stack
    stack_len = length(stack_stack)
    @inbounds for nt = left:left + len - 1
        count = stack_stack[genome_wrap(stack_len, nt)]
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
        depth = sum_count / (max_count * len)
    end
    return depth, coverage / len
end

function alignTemplateToStack(feature_stack::FeatureStack, shadowstack::ShadowStack)::Tuple{Int32,Int32}
    stack = feature_stack.stack
    glen = length(stack)
    median_length = feature_stack.template.median_length
    score = 0
    @inbounds for nt = 1:median_length
        gw = genome_wrap(glen, nt)
        score += stack[gw] + shadowstack[gw]
    end
    max_score = score
    best_hit = 1

    @inbounds for nt = 2:glen
        st = shadowstack[nt]
        score -= stack[nt - 1] + st
        score += stack[genome_wrap(glen, nt + median_length - 1)] + st
        if score > max_score
            max_score = score
            best_hit = nt
        end
    end
    max_score <= 0 && return 0, 0
    return best_hit, median_length
end

function getFeaturePhaseFromAnnotationOffsets(feat::Feature, annotations::Vector{Annotation})::Int8
    phases = Int8[]
    for annotation in annotations[findall(x -> x.path == feat.path, annotations)]
        if rangesOverlap(feat.start, feat.length, annotation.start, annotation.length)
            # estimating feature phase from annotation phase
            phase = phaseCounter(annotation.phase, feat.start - annotation.start)
            push!(phases, phase)
        end
    end
    if length(phases) == 0
        @warn "No annotations found for $(feat.path)"
        return zero(Int8)
    end
    return StatsBase.mode(phases) # return most common phase
end


function weightedMode(values::Vector{Int32}, weights::Vector{Float32})::Int32
    w, v = findmax(StatsBase.addcounts!(Dict{Int32,Float32}(), values, weights))
    return v
end

# uses weighted mode, weighting by alignment length and distance from boundary
function refineMatchBoundariesByOffsets!(feat::Feature, annotations::Vector{Annotation}, 
            target_length::Integer, coverages::Dict{String,Float32})  # Tuple{Feature,Vector{Int32},Vector{Float32}}
    # grab all the matching features
    matching_annotations = findall(x -> x.path == feat.path, annotations)
    isempty(matching_annotations) && return #  feat, [], []
    overlapping_annotations = Annotation[]
    minstart = target_length
    maxend = 1
    for annotation in annotations[matching_annotations]
        annotation.offset5 >= feat.length && continue
        annotation.offset3 >= feat.length && continue
        if rangesOverlap(feat.start, feat.length, annotation.start, annotation.length)
            minstart = min(minstart, annotation.start)
            maxend = max(maxend, annotation.start + annotation.length - 1)
            push!(overlapping_annotations, annotation)
        end
    end
    n = length(overlapping_annotations)
    end5s = Vector{Int32}(undef, n)
    end5ws = Vector{Float32}(undef, n)
    end3s = Vector{Int32}(undef, n)
    end3ws = Vector{Float32}(undef, n)

    @inbounds for (i, annotation) in enumerate(overlapping_annotations)
        # predicted 5' end is annotation start - offset5
        end5s[i] = annotation.start - annotation.offset5
        # weights
        coverage = coverages[annotation.genome_id]
        weight = ((feat.length - (annotation.start - minstart)) / feat.length) * coverage
        end5ws[i] =  weight * weight
        # predicted 3' end is feature start + feature length + offset3 - 1
        end3s[i] = annotation.start + annotation.length + annotation.offset3 - 1
        # weights
        weight = ((feat.length - (maxend - (annotation.start + annotation.length - 1))) / feat.length) * coverage
        end3ws[i] = weight * weight
    end
    left = feat.start
    right = feat.start + feat.length - 1
    # println(feat.path)
    if !isempty(end5s)
        left = weightedMode(end5s, end5ws)
    end
    if !isempty(end3s)
        right = weightedMode(end3s, end3ws)
    end
    feat.start = genome_wrap(target_length, left)
    feat.length = right - left + 1
    # return feat # , end5s, end5ws
end

function groupFeaturesIntoGeneModels(features::AFeature)::AAFeature
    gene_models = AAFeature()
    current_model = Feature[]
    for feature in features
        if isempty(current_model)
            push!(current_model, feature)
        elseif getFeatureName(current_model[1]) == getFeatureName(feature)
            push!(current_model, feature)
        else
            sort!(current_model, by=x -> x.start)
            push!(gene_models, current_model)
            current_model = Feature[]
            push!(current_model, feature)
        end
    end
    if length(current_model) > 0
        sort!(current_model, by=x -> x.start)
        push!(gene_models, current_model)
    end
    return gene_models
end

function translateFeature(genome::DNAString, feat::Feature)
    @assert feat.start > 0
    feat.length < 3 && return ""

    fend = feat.start + feat.length - 3
    if fend > length(genome) - 2
        fend = length(genome) - 2
        @warn "translateFeature: feature points past genome"
    end
    translateDNA(genome, feat.start + feat.phase, fend)
end

function translateModel(genomeloop::DNAString, model::AFeature)::String
    DNA = Vector{String}()
    sizehint!(DNA, length(model))
    for (i, feat) in enumerate(model)
        getFeatureType(feat) ≠ "CDS" && continue
        start = feat.start
        if i == 1
            start += feat.phase
        end
        push!(DNA, SubString(genomeloop, start, feat.start + feat.length - 1))
    end
    
    dna = join(DNA, "")
    translateDNA(dna)
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

function findStartCodon2!(cds::Feature, genome_length::Integer, genomeloop::DNAString, predicted_starts,
                          predicted_starts_weights, feature_stack, shadowstack)
    # define range in which to search from upstream stop to end of feat
    start5 = genome_wrap(genome_length, cds.start + cds.phase)
    codon = SubString(genomeloop, start5, start5 + 2)
    while !isStopCodon(codon, false)
        start5 = genome_wrap(genome_length, start5 - 3)
        codon = SubString(genomeloop, start5, start5 + 2)
    end
    search_range = start5:cds.start + cds.length - 1
    window_size = length(search_range)

    codons = zeros(window_size)
    phase = zeros(window_size)
    cumulative_stack_coverage = zeros(window_size)
    cumulative_shadow_coverage = zeros(window_size)

    stack_coverage = 0
    shadow_coverage = 0
    for (i, nt) in enumerate(search_range)
        gw = genome_wrap(genome_length, nt)
        codon = SubString(genomeloop, nt, nt + 2)
        codons[i] = startScore(cds, codon)
        if i % 3 == 1
            phase[i] = 1.0
        end
        if feature_stack[gw] > 0
            stack_coverage += 1
        elseif shadowstack[gw] < 0
            shadow_coverage -= 1
        end
        if stack_coverage > 3
            cumulative_stack_coverage[i] = 3 - stack_coverage
        end
        cumulative_shadow_coverage[i] = shadow_coverage
    end

    predicted_starts = zeros(window_size)
    for (s, w) in zip(predicted_starts, predicted_starts_weights)
        for (i, nt) in enumerate(search_range)
            distance = 1 + (s - nt)^2
            predicted_starts[i] = w / distance
        end
    end

    # combine vectors using parameter weights
    result = codons .* phase .* cumulative_stack_coverage .+ cumulative_shadow_coverage .+ predicted_starts
    # set feat.start to highest scoring position
    maxstartscore, maxstartpos = findmax(result)
    maxstartpos += start5 - 1
    cds.length += cds.start - maxstartpos
    cds.start = genome_wrap(genome_length, maxstartpos)
    # set phase to zero
    cds.phase = 0
    return cds
end

function findStartCodon!(cds::Feature, genome_length::Int32, genomeloop::DNAString)
    # assumes phase has been correctly set
    # search for start codon 5'-3' beginning at cds.start, save result; abort if stop encountered
    start3 = cds.start + cds.phase
    codon = SubString(genomeloop, start3, start3 + 2)
    if isStartCodon(codon, true, true)  # allow ACG and GTG codons if this is predicted start
        cds.start = start3
        cds.length -= cds.phase
        cds.phase = 0;
        return cds
    end

    allowACG = false
    allowGTG = false
    if getFeatureName(cds) == "rps19"
        allowGTG = true
    end
    
    start3 = nothing
    for i in iter_wrap(cds.start + cds.phase:3:cds.start + cds.phase + genome_length - 3, genome_length)
        codon = SubString(genomeloop, i, i + 2)
        if isStartCodon(codon, allowACG, allowGTG)
            start3 = i >= cds.start + cds.phase ? i : i + genome_length
            break
        end
        if isStopCodon(codon, false)
            break
        end
    end
    # search for start codon 3'-5' beginning at cds.start, save result; abort if stop encountered
    start5 = nothing
    for i in iter_wrap(cds.start + cds.phase:-3:cds.start + cds.phase - genome_length + 3, genome_length)
        codon = SubString(genomeloop, i, i + 2)
        if isStartCodon(codon, allowACG, allowGTG)
            start5 = i <= cds.start + cds.phase ? i : i - genome_length
            break
        end
        if isStopCodon(codon, false)
            break
        end
    end
    # return cds with start set to nearer of the two choices
    if isnothing(start3) && isnothing(start5)
        @warn "Couldn't find start for $(cds.path)"
        return cds
    end
    if isnothing(start3)
        cds.length += cds.start - start5
        cds.start = genome_wrap(genome_length, start5)
    elseif isnothing(start5)
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

function setLongestORF!(feat::Feature, genome_length::Int32, targetloop::DNAString)
    orfs = []
    translation_start = feat.start
    translation_stop = translation_start + 2
    fend = feat.start + feat.length - 3

    for nt = translation_start + feat.phase:3:fend
        translation_stop = nt + 2
        codon = SubString(targetloop, nt, translation_stop)
        if isStopCodon(codon, false)
            push!(orfs, (translation_start, translation_stop, true))
            translation_start = nt + 3
        end
    end
    push!(orfs, (translation_start, translation_stop, false))
    maxgap = 0
    longest_orf = nothing
    for orf in orfs
        gap = orf[2] - orf[1] + 1
        if gap > maxgap
            maxgap = gap
            longest_orf = orf
        end
    end
    longest_orf === nothing && return

    if longest_orf[1] > feat.start # must be an internal stop
        feat.phase = 0
    end
    feat.length = longest_orf[2] - longest_orf[1] + 1
    feat.start = genome_wrap(genome_length, longest_orf[1])

    if !longest_orf[3]
        # orf is still open, so extend until stop
        translation_start = longest_orf[2] - 2
        codon = SubString(targetloop, translation_start, translation_start + 2)
        if isStopCodon(codon, true)
            return feat
        end
        for nt in iter_wrap(translation_start:3:translation_start + genome_length - 3, genome_length)
        # for nt = translation_start:3:length(targetloop) - 3 # sic! - 3
            codon = SubString(targetloop, nt, nt + 2)
            if isStopCodon(codon, false)
                break
            end
            feat.length += 3
        end
    end
end

function refineBoundariesbyScore!(feat1::Feature, feat2::Feature, genome_length::Int32, stacks::DFeatureStack)
    # feat1 should be before feat2
    start = feat1.start + feat1.length - 1
    range_to_test = min(start, feat2.start):max(start, feat2.start)
    
    feat1_stack = stacks[feat1.path].stack
    feat2_stack = stacks[feat2.path].stack

    # score = sum(@view feat2_stack[range_to_test]) - sum(@view feat1_stack[range_to_test])
    score = sum(iter_wrap(range_to_test, feat2_stack)) - sum(iter_wrap(range_to_test, feat1_stack))
    maxscore = score
    fulcrum = first(range_to_test) - 1
    @inbounds for i in range_to_test
        gw = genome_wrap(genome_length, i)
        score += 2 * feat1_stack[gw]
        score -= 2 * feat2_stack[gw]
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

function refineGeneModels!(gene_models::AAFeature, genome_length::Int32, targetloop::DNAString,
                          annotations::Vector{Annotation},
                          feature_stacks::DFeatureStack)::AAFeature
    last_cds_examined = nothing
    for model in gene_models
        isempty(model) && continue
        # sort features in model by mid-point to avoid cases where long intron overlaps short exon
        sort!(model, by=f -> f.start + f.length / 2)
        last_exon = last(model)
        # if CDS, find phase, stop codon and set feature.length
        if isType(last_exon, "CDS")
            last_exon.phase = getFeaturePhaseFromAnnotationOffsets(last_exon, annotations)

            if !isFeatureName(last_exon, "rps12A")
                setLongestORF!(last_exon, genome_length, targetloop)
            end
            last_cds_examined = last_exon
        end
        
        @inbounds for i in length(model) - 1:-1:1
            feature = model[i]
            # check adjacent to last exon, if not...
            gap = model[i + 1].start - (feature.start + feature.length)
            if gap ≠ 0 && gap < 100
                refineBoundariesbyScore!(feature, model[i + 1], genome_length, feature_stacks)
                gap = model[i + 1].start - (feature.start + feature.length)
            end
            if gap ≠ 0
                @warn "Non-adjacent boundaries for $(feature.path) $(model[i + 1].path)"
            end
            # if CDS, check phase is compatible
            if isType(feature, "CDS")
                feature.phase = getFeaturePhaseFromAnnotationOffsets(feature, annotations)
                if last_cds_examined !== nothing && phaseCounter(feature.phase, feature.length % Int32(3)) != last_cds_examined.phase
                    @warn "Incompatible phases for $(feature.path) $(last_cds_examined.path)"
                    # refine boundaries (how?)
                end
                last_cds_examined = feature
            end
        end
        # if CDS, find start codon and set feature.start
        first_exon = first(model)
        if isType(first_exon, "CDS")
            first_exon.phase = getFeaturePhaseFromAnnotationOffsets(first_exon, annotations)
            if !isFeatureName(last_exon, "rps12B")
                findStartCodon!(first_exon, genome_length, targetloop)
                # first_exon = findStartCodon2!(first_exon,genome_length,targetloop)
            end
        end
    end
    return gene_models
end



struct SFF
    gene::String
    feature::Feature
    gene_length::Int
    exon_count::Int
    avg::Float64
    depth::Float64
    coverage::Float64
    hasStart::Bool
    hasPrematureStop::Bool
end

gene_length(model::Vector{Feature}) = begin
    maximum([f.start + f.length for f in model]) - minimum([f.start for f in model])
end

gene_length(model::Vector{SFF}) = begin
    maximum([m.feature.start + m.feature.length for m in model]) - minimum([m.feature.start for m in model])
end

function toSFF(model::Vector{Feature}, targetloop::DNAString, feature_stacks::DFeatureStack)::Vector{SFF}
    first_model = first(model)
    gene = getFeatureName(first_model)
    ret = []
    exon_count = 0
    cds = false
    for f in model
        t = getFeatureType(f)
        if t == "CDS"
            cds = true
        end
        if t ≠ "intron"
            exon_count += 1
        end
    end
    hasStart = true
    hasPrematureStop = false
    if cds
        protein = translateModel(targetloop, model)
        if gene ≠ "rps12B" && !isStartCodon(SubString(targetloop, first_model.start, first_model.start + 2), true, true)
            hasStart = false
        end
        stop_position = findfirst(isequal('*'), protein)
        if !isnothing(stop_position) && stop_position < length(protein)
            hasPrematureStop = true
        end
    end
    genelength = gene_length(model)
    genelength <= 0 && return ret
    
    for feature in model
        # stack = feature_stacks[findfirst(x -> x.path == f.path, feature_stacks)]
        stack = get(feature_stacks, feature.path, nothing)
        if stack === nothing
            @error "$(gene): no stack for path $(f.path)"
            continue
        end

        depth, coverage = getDepthAndCoverage(stack, feature.start, feature.length)
        depth == 0 && continue
        coverage == 0 && continue
        avg = feature.length / stack.template.median_length
        push!(ret, SFF(gene, feature, genelength, exon_count, avg, depth, coverage, hasStart, hasPrematureStop))
    end
    return ret
end

# ChloeIO = Union{IOStream,IOBuffer,GZipStream}

function writeModelToSFF(outfile::IO, model_id::String,
                        sffs::Vector{SFF},
                        gene_exons::Dict{String,Int32},
                        maxlengths::Dict{String,Int32},
                        strand::Char)
    
    for sff in sffs
        # stack = feature_stacks[findfirst(x -> x.path == f.path, feature_stacks)]
        feature = sff.feature
        write(outfile, model_id)
        write(outfile, "/")
        write(outfile, getFeatureSuffix(feature))
        write(outfile, "\t")
        write(outfile, join([strand,string(feature.start),string(feature.length),string(feature.phase)], "\t"))
        write(outfile, "\t")
        write(outfile, join([@sprintf("%.3f",sff.avg),@sprintf("%.3f",sff.depth),@sprintf("%.3f",sff.coverage)], "\t"))
        write(outfile, "\t")

        pseudo = "possible pseudogene"
        if sff.gene_length - get(maxlengths, sff.gene, 0) < -10
            write(outfile, "$(pseudo), shorter than 2nd copy")
            pseudo = ""
        end
        expected_exons = get(gene_exons, sff.gene, 0)
        if sff.exon_count < expected_exons
            write(outfile, "$(pseudo), missing exon(s)")
            pseudo = ""
        end
        if !sff.hasStart
            write(outfile, "$(pseudo), no start codon")
            pseudo = ""
        end
        if sff.hasPrematureStop
            write(outfile, "$(pseudo), premature stop codon")
        end

        # if #too short compared to template
        #     write(outfile,"possible pseudogene, too short ")
        # end
        # if #too short divergent compared to template
        #     write(outfile,"possible pseudogene, too divergent ")
        # end
        write(outfile, "\n")
    end
end

function calc_maxlengths(models::FwdRev{Vector{Vector{SFF}}})::Dict{String,Int32}

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

MaybeIR = Union{AlignedBlock,Nothing}

function writeSFF(outfile::Union{String,IO},
                id::String, # NCBI id
                genome_length::Int32,
                gene_exons::Dict{String,Int32}, # reference exons lengths
                models::FwdRev{Vector{Vector{SFF}}},
                ir::MaybeIR=nothing)


    function getModelID!(model_ids::Dict{String,Int32}, model::Vector{SFF})
        gene_name = first(model).gene
        instance_count::Int32 = 1
        model_id = "$(gene_name)/$(instance_count)"
        while get(model_ids, model_id, 0) ≠ 0
            instance_count += 1
            model_id = "$(gene_name)/$(instance_count)"
        end
        model_ids[model_id] = instance_count
        return model_id
    end
    

    maxlengths = calc_maxlengths(models)


    function out(outfile::IO)
        model_ids = Dict{String,Int32}()
        write(outfile, id, "\t", string(genome_length), "\n")
        for model in models.forward
            isempty(model) && continue
            model_id = getModelID!(model_ids, model)
            writeModelToSFF(outfile, model_id, model, gene_exons, maxlengths, '+')
        end
        for model in models.reverse
            isempty(model) && continue
            model_id = getModelID!(model_ids, model)
            writeModelToSFF(outfile, model_id, model, gene_exons, maxlengths, '-')
        end
        if ir !== nothing
            println(outfile, "IR/1/repeat_region/1\t+\t$(ir[1])\t$(ir[3])\t0\t0\t0\t0\t")
            println(outfile, "IR/2/repeat_region/1\t-\t$(ir[2])\t$(ir[3])\t0\t0\t0\t0\t")
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
