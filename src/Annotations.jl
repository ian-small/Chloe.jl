using Printf

mutable struct Feature
    path::String
    start::Int32
    length::Int32
    # phase is the number of nucleotides to skip at the start of the sequence 
    # to be in the correct reading frame
    phase::Int8
    _path_components::Array{String}
    Feature(path, start, length, phase) = new(path, start, length, phase, split(path, '/'))
end
const featurePath(feature::Feature) = feature._path_components # split(feat.path, '/')
function getFeatureName(feature::Feature)
    feature._path_components[1]
end
function getFeatureType(feature::Feature)
    feature._path_components[3]
end

function isType(feature::Feature, gene::String)
    gene == getFeatureType(feature)
end
function isFeatureName(feature::Feature, name::String)
    name == getFeatureName(feature)
end

AFeature = Array{Feature}
AAFeature = Array{AFeature}
# entire set of Features for one strand of one genome
struct FeatureArray
    genome_id::String
    genome_length::Int32
    strand::Char
    features::AFeature
end

function readFeatures(file::String)::Tuple{FeatureArray,FeatureArray}
    open(file) do f
        header = split(readline(f), '\t')
        f_strand_features = FeatureArray(header[1], parse(Int32, header[2]), '+', AFeature(undef, 0))
        r_strand_features = FeatureArray(header[1], parse(Int32, header[2]), '-', AFeature(undef, 0))
        while !eof(f)
            fields = split(readline(f), '\t')
            feature = Feature(fields[1], parse(Int, fields[3]), parse(Int, fields[4]), parse(Int, fields[5]))
            if fields[2][1] == '+'
                push!(f_strand_features.features, feature)
            else
                push!(r_strand_features.features, feature)
            end
        end
        return f_strand_features, r_strand_features
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

# checks all blocks for overlap so could be speeded up by using an interval tree
function addOverlapBlocks(genome_id::String, feature::Feature, blocks::AlignedBlocks)::Array{Annotation}
    pushed_features = Array{Annotation}(undef, 0)
    feature_type = getFeatureType(feature)
    for block in blocks
        if rangesOverlap(feature.start, feature.length, block[1], block[3])
            startA = max(feature.start, block[1])
            flength = min(feature.start + feature.length, block[1] + block[3]) - startA
            flength <= 0 && continue
            offset5 = startA - feature.start
            if (feature_type == "CDS")
                phase = phaseCounter(feature.phase, offset5)
            else
                phase = 0
            end
            path_components = featurePath(feature)
            annotation_path = join([path_components[1],"?",path_components[3],path_components[4]], "/")
            pushed_feature = Annotation(genome_id, annotation_path, 
                                        startA - block[1] + block[2],
                                        flength, offset5, feature.length - (startA - feature.start) - flength, phase)

            push!(pushed_features, pushed_feature)
        end
    end
    return pushed_features
end

struct FeatureTemplate
    path::String  # similar to .sff path
    threshold_counts::Float32
    threshold_coverage::Float32
    median_length::Int32
end

using StatsBase

function readTemplates(file::String)::Tuple{Array{FeatureTemplate},Dict{String,Int32}}
    templates = FeatureTemplate[]
    gene_exons = String[]
    if !isfile(file)
        error("$(file) is not a file")
    end

    open(file) do f
        header = readline(f)
        while !eof(f)
            fields = split(readline(f), '\t')
            path_components = split(fields[1], '/')
            if path_components[3] ≠ "intron"
                push!(gene_exons, path_components[1])
            end
            template = FeatureTemplate(fields[1], 
                parse(Float32, fields[2]),
                parse(Float32, fields[3]),
                parse(Int32, fields[4]))
            push!(templates, template)
        end
    end
    return sort!(templates, by = x->x.path), countmap(gene_exons)
end

struct AnnotationArray
    genome_id::String
    strand::Char
    annotations::Array{Annotation}
end

function findOverlaps(ref_featurearray::FeatureArray, aligned_blocks::AlignedBlocks)::Array{Annotation}
    annotations = Array{Annotation}(undef, 0)
    for feature in ref_featurearray.features
        new_annotations = addOverlapBlocks(ref_featurearray.genome_id, feature, aligned_blocks)
        if !isempty(new_annotations)
            # pushed_features.annotations = cat(pushed_features.annotations, new_features, dims = 1)
            annotations = cat(annotations, new_annotations, dims = 1)
        end
    end
    return annotations
end

struct FeatureStack
    path::String
    stack::Array{Int32}
    template::FeatureTemplate
end

AFeatureStack = Array{FeatureStack}
ShadowStack = Array{Int32}

function fillFeatureStack(target_length::Int32, annotations::AnnotationArray,
    templates::Array{FeatureTemplate})::Tuple{AFeatureStack,ShadowStack}
    stacks = AFeatureStack(undef, 0)
    shadowstack::ShadowStack = fill(Int32(-1), target_length) # will be negative image of all stacks combined,
    # initialised to small negative number; acts as prior expectation for feature-finding
    for annotation in annotations.annotations
        template_index = findfirst(x->x.path == annotation.path, templates)
        if template_index === nothing
            @error "Can't find template for $(annotation.path)"
        end
        if isempty(stacks) || (index = findfirst(x->x.path == annotation.path, stacks)) === nothing
            if template_index === nothing
                continue # what to do?
            end
            stack = FeatureStack(annotation.path, zeros(Int32, target_length), templates[template_index])
            push!(stacks, stack)
        else
            stack = stacks[index]
        end
        for i = annotation.start:annotation.start + annotation.length - Int32(1)
            gw = genome_wrap(target_length, i)
            stack.stack[gw] += 3 # +1 to counteract shadowstack initialisation, +1 to counteract addition to shadowstack
            shadowstack[gw] -= 1
        end
    end
    @debug "found $(length(stacks)) FeatureStacks from $(length(annotations.annotations)) annotations"
    return stacks, shadowstack
end

function expandBoundaryInChunks(feature_stack::FeatureStack, shadowstack::ShadowStack, origin, direction, max)::Integer
    glen = length(shadowstack)
    pointer = origin

    # expand in chunks
    for chunksize in [100,80,60,40,30,20,10,5,1]
        for i = direction:direction * chunksize:direction * max
            sumscore = 0
            for j = 1:chunksize
                sumscore += feature_stack.stack[genome_wrap(glen, pointer + direction * j)] + shadowstack[genome_wrap(glen, pointer + direction * j)]
            end
            sumscore < 0 && break
            pointer += direction * chunksize
        end
    end
    return genome_wrap(glen, pointer)
end

function expandBoundary(feature_stack::FeatureStack, shadowstack::ShadowStack, origin, direction, max)
    glen = length(shadowstack)
    pointer = origin
    for i = direction:direction:direction * max
        sumscore = 0
        sumscore += feature_stack.stack[genome_wrap(glen, pointer + direction)] + shadowstack[genome_wrap(glen, pointer + direction)]
        sumscore < 0 && break
        pointer += direction
    end
    return pointer # not wrapped
end

function getDepthAndCoverage(feature_stack::FeatureStack, left::Int32, len::Int32)::Tuple{Float64,Float64}
    coverage = 0
    max_count = 0
    sum_count = 0
    for nt = left:left + len - 1
        count = feature_stack.stack[genome_wrap(length(feature_stack.stack), nt)]
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
    tlen = feature_stack.template.median_length
    score = 0
    best_hit = 1
    for nt = 1:tlen
        count = stack[genome_wrap(glen, nt)] + shadowstack[genome_wrap(glen, nt)]
        score += count
    end
    max_score = score
    for nt = 2:glen
        count = stack[genome_wrap(glen, nt - 1)] + shadowstack[genome_wrap(glen, nt)]
        score -= count
        count = stack[genome_wrap(glen, nt + tlen - 1)] + shadowstack[genome_wrap(glen, nt)]
        score += count
        if score > max_score
            max_score = score
            best_hit = nt
        end
    end
    max_score <= 0 && return 0, 0
    return best_hit, tlen
end

function getModelID!(model_ids::Dict{String,Int32}, model::AFeature)
    gene_name = getFeatureName(first(model))
    instance_count = 1
    model_id = gene_name * "/" * string(instance_count)
    while get(model_ids, model_id, 0) ≠ 0
        instance_count += 1
        model_id = gene_name * "/" * string(instance_count)
    end
    push!(model_ids, model_id => instance_count)
    return model_id
end

function writeSFF(outfile::String, fstrand_features::FeatureArray, rstrand_features::FeatureArray)
    open(outfile, "w") do outfile
        write(outfile, fstrand_features.genome_id, "\t", string(fstrand_features.genome_length), "\n")
        for f in fstrand_features.features
            write(outfile, f.path)
            write(outfile, "\t")
            write(outfile, fstrand_features.strand)
            write(outfile, "\t")
            write(outfile, string(f.start))
            write(outfile, "\t")
            write(outfile, string(f.length))
            write(outfile, "\n")
        end
        for f in rstrand_features.features
            write(outfile, f.path)
            write(outfile, "\t")
            write(outfile, rstrand_features.strand)
            write(outfile, "\t")
            write(outfile, string(f.start))
            write(outfile, "\t")
            write(outfile, string(f.length))
            write(outfile, "\n")
        end
    end
end

function getFeaturePhaseFromAnnotationOffsets(feat::Feature, annotations::AnnotationArray)::Int8
    matching_annotations = findall(x->x.path == feat.path, annotations.annotations)
    phases = Int8[]
    for annotation in annotations.annotations[matching_annotations]
        if rangesOverlap(feat.start, feat.length, annotation.start, annotation.length)
            # estimating feature phase from annotation phase
            phase = phaseCounter(annotation.phase, feat.start - annotation.start)
            push!(phases, phase)
        end
    end
    if length(phases) == 0
        @warn "No annotations found for $(feat.path)"
        return 0
    end
    return StatsBase.mode(phases) # return most common phase
end

function weightedMode(values::Array{Int32}, weights::Array{Float32})::Int32
    weightedCounts = zeros(0, 2)
    for (v, w) in zip(values, weights)
        row = findfirst(isequal(v), weightedCounts[:,1])
        if isnothing(row)
            weightedCounts = vcat(weightedCounts, [v w])
        else
            weightedCounts[row,2] += w
        end
    end
    maxweight, maxrow = findmax(weightedCounts[:,2])
    return Int32(weightedCounts[maxrow,1])
end

# uses weighted mode, weighting by alignment length and distance from boundary
function refineMatchBoundariesByOffsets!(feat::Feature, annotations::AnnotationArray, 
            target_length::Integer, coverages::Dict{String,Float32})::Tuple{Feature,Array{Int32},Array{Float32}}
    # grab all the matching features
    matching_annotations = findall(x->x.path == feat.path, annotations.annotations)
    isempty(matching_annotations) && return feat, [], []
    overlapping_annotations = []
    minstart = target_length
    maxend = 1
    for annotation in annotations.annotations[matching_annotations]
        annotation.offset5 >= feat.length && continue
        annotation.offset3 >= feat.length && continue
        if rangesOverlap(feat.start, feat.length, annotation.start, annotation.length)
            if annotation.start < minstart; minstart = annotation.start; end
            if annotation.start + annotation.length - 1 > maxend; maxend = annotation.start + annotation.length - 1; end
            push!(overlapping_annotations, annotation)
        end
    end
    end5s = Array{Int32}(undef, 0)
    end5ws = Array{Float32}(undef, 0)
    end3s = Array{Int32}(undef, 0)
    end3ws = Array{Float32}(undef, 0)

    for annotation in overlapping_annotations
        # predicted 5' end is annotation start - offset5
        push!(end5s, annotation.start - annotation.offset5)
        # weights
        coverage = coverages[annotation.genome_id]
        weight = ((feat.length - (annotation.start - minstart)) / feat.length) * coverage
        push!(end5ws, weight * weight)
        # predicted 3' end is feature start + feature length + offset3 - 1
        push!(end3s, annotation.start + annotation.length + annotation.offset3 - 1)
        # weights
        weight = ((feat.length - (maxend - (annotation.start + annotation.length - 1))) / feat.length) * coverage
        push!(end3ws, weight * weight)
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
    return feat, end5s, end5ws
end




function groupFeaturesIntoGeneModels(features::FeatureArray)::AAFeature
    gene_models = AAFeature(undef, 0)
    current_model = Feature[]
    for feature in features.features
        if isempty(current_model)
            push!(current_model, feature)
        elseif getFeatureName(current_model[1]) == getFeatureName(feature)
            push!(current_model, feature)
        else
            sort!(current_model, by = x->x.start)
            push!(gene_models, current_model)
            current_model = []
            push!(current_model, feature)
        end
    end
    sort!(current_model, by = x->x.start)
    push!(gene_models, current_model)
    return gene_models
end

function translateFeature(genome::DNAString, feat::Feature)
    @assert feat.start > 0
    feat.length < 3 && return ""

    peptide = Array{Char}(undef, fld(feat.length, 3))

    aa = 0

    for i = (feat.start + feat.phase):3:feat.start + feat.length - 3
        aa += 1
        peptide[aa] = get(genetic_code, SubString(genome, i, i + 2), 'X')
    end
    return String(peptide)
end

function translateModel(genome::DNAString, model::AFeature)

    DNA = ""
    for (i, feat) in enumerate(model)
        getFeatureType(feat) ≠ "CDS" && continue
        start = feat.start
        if i == 1
            start += feat.phase
        end
        DNA = DNA * SubString(genome, start, feat.start + feat.length - 1)
    end

    peptide = Array{Char}(undef, fld(length(DNA), 3))

    aa = 0
    for i = 1:3:length(DNA) - 2
        aa += 1
        peptide[aa] = get(genetic_code, SubString(DNA, i, i + 2), 'X')
    end
    return String(peptide)
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
    start5 = cds.start + cds.phase
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
    i = 1
    stack_coverage = 0
    shadow_coverage = 0
    for nt in search_range
        codon = SubString(genomeloop, nt, nt + 2)
        codons[i] = startScore(codon)
        if i % 3 == 1; phase[i] = 1.0; end
        if feature_stack[nt] > 0
            stack_coverage += 1
        elseif shadowstack[nt] < 0
            shadow_coverage -= 1
        end
        if stack_coverage > 3; cumulative_stack_coverage[i] = 3 - stack_coverage; end
        cumulative_shadow_coverage[i] = shadow_coverage
        i += 1
    end

    predicted_starts = zeros(window_size)
    for (s, w) in zip(predicted_starts, predicted_starts_weights)
        i = 1
        for nt in search_range
            distance = 1 + (s - nt)^2
            predicted_starts[i] = w / distance
            i += 1
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
    Three = Int32(3)
    start3 = cds.start + cds.phase
    codon = SubString(genomeloop, start3, start3 + 2)
    if isStartCodon(codon, true, true)  # allow ACG and GTG codons if this is predicted start
        cds.start = start3
        cds.length -= cds.phase
        cds.phase = 0;
        return cds
    end

    gene = getFeatureName(cds)
    allowACG = false
    allowGTG = false
    if gene == "rps19"; allowGTG = true; end

    while (!isStartCodon(codon, allowACG, allowGTG) && start3 <= genome_length)
        if isStopCodon(codon, false)
            start3 = nothing
            break
        end
        start3 = genome_wrap(genome_length, start3 + Three)
        codon = SubString(genomeloop, start3, start3 + 2)
    end
    # search for start codon 3'-5' beginning at cds.start, save result; abort if stop encountered
    start5 = cds.start + cds.phase
    codon = SubString(genomeloop, start5, start5 + 2)
    while (!isStartCodon(codon, allowACG, allowGTG) && start5 > 0)
        if isStopCodon(codon, false)
            start5 = nothing
            break
        end
        start5 = genome_wrap(genome_length, start5 - Three)
        codon = SubString(genomeloop, start5, start5 + 2)
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
    zero, one, two, three = Int32(0), Int32(1), Int32(2), Int32(3)
    translation_start = feat.start
    translation_stop = translation_start + two
    for nt = translation_start + feat.phase:three:feat.start + feat.length - three
        translation_stop = nt + two
        codon = SubString(targetloop, nt, translation_stop)
        if isStopCodon(codon, false)
            push!(orfs, (translation_start, translation_stop, true))
            translation_start = nt + three
        end
    end
    push!(orfs, (translation_start, translation_stop, false))
    maxgap = zero
    local longest_orf
    for orf in orfs
        gap = orf[2] - orf[1] + one
        if gap > maxgap
            maxgap = gap
            longest_orf = orf
        end
    end

    if longest_orf[1] > feat.start # must be an internal stop
        feat.phase = 0
    end
    feat.length = longest_orf[2] - longest_orf[1] + one
    feat.start = genome_wrap(genome_length, longest_orf[1])

    if !longest_orf[3]
       # orf is still open, so extend until stop
        translation_start = longest_orf[2] - two
        codon = SubString(targetloop, translation_start, translation_start + two)
        if isStopCodon(codon, true)
            return feat
        end
        for nt = translation_start:3:length(targetloop) - 3
            codon = SubString(targetloop, nt, nt + 2)
            if isStopCodon(codon, false)
                break
            end
            feat.length += 3
        end
    end
end

function refineBoundariesbyScore!(feat1::Feature, feat2::Feature, stacks::Array{FeatureStack})
    # feat1 shoud be before feat2
    range_to_test = min(feat1.start + feat1.length - 1, feat2.start):max(feat1.start + feat1.length - 1, feat2.start)

    feat1_stack = stacks[findfirst(x->x.path == feat1.path, stacks)]
    feat2_stack = stacks[findfirst(x->x.path == feat2.path, stacks)]
    score = sum(feat2_stack.stack[range_to_test]) - sum(feat1_stack.stack[range_to_test])
    maxscore = score
    fulcrum = first(range_to_test) - 1
    for i in range_to_test
        score += 2 * feat1_stack.stack[genome_wrap(length(feat1_stack.stack), i)]
        score -= 2 * feat2_stack.stack[genome_wrap(length(feat2_stack.stack), i)]
        if score > maxscore
            maxscore = score
            fulcrum = i
        end
    end
    # fulcrum should point to end of feat1
    feat1.length = fulcrum - feat1.start + 1
    # fulcrum+1 should point to start of feat2
    feat2.length += feat2.start - (fulcrum + 1)
    feat2.start = fulcrum + 1
end

function refineGeneModels!(gene_models::AAFeature, genome_length::Int32, targetloop::DNAString,
                          annotations::AnnotationArray,
                          feature_stacks::AFeatureStack)::AAFeature
    for model in gene_models
        isempty(model) && continue
        # sort features in model by mid-point to avoid cases where long intron overlaps short exon
        sort!(model, by = f->f.start + f.length / 2)
        # @debug "model"  model
        last_exon = last(model)
        # if CDS, find phase, stop codon and set feature.length
        if isType(last_exon, "CDS")
            translation = translateFeature(targetloop, last_exon)
            # @debug "translation" translation
            last_exon.phase = getFeaturePhaseFromAnnotationOffsets(last_exon, annotations)

            if !isFeatureName(last_exon, "rps12A")
                setLongestORF!(last_exon, genome_length, targetloop)
            end

            translation = translateFeature(targetloop, last_exon)

            last_cds_examined = last_exon
        end
        for i in length(model) - 1:-1:1
            feature = model[i]
            # check adjacent to last exon, if not...
            gap = model[i + 1].start - (feature.start + feature.length)
            if gap ≠ 0 && gap < 100
                refineBoundariesbyScore!(feature, model[i + 1], feature_stacks)
                gap = model[i + 1].start - (feature.start + feature.length)
            end
            if gap ≠ 0
                @warn "Non-adjacent boundaries for $(feature.path) $(model[i + 1].path)"
            end
            # @debug "feature" feature
            # if CDS, check phase is compatible
            if isType(feature, "CDS")
                feature.phase = getFeaturePhaseFromAnnotationOffsets(feature, annotations)
                if (@isdefined last_cds_examined) && phaseCounter(feature.phase, feature.length % Int32(3)) != last_cds_examined.phase
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
                first_exon = findStartCodon!(first_exon, genome_length, targetloop)
                # first_exon = findStartCodon2!(first_exon,genome_length,targetloop)
            end
        end
    end
    return gene_models
end

MaybeFeature = Union{Feature,Nothing}
MaybeAFeature = Union{AFeature,Nothing}
function getFeatureByName(fname::String, features::FeatureArray)::MaybeFeature
    for feat in features.features
        if isFeatureName(feat, fname)
            return feat
        end
    end
    return nothing
end

function getGeneModelByName(gm_name::String, gene_models::AAFeature)::MaybeAFeature
    for model in gene_models
        if isempty(model)
            continue
        end
        if isFeatureName(model[1], gm_name)
            return model::AFeature
        end
    end
    return nothing
end

function writeModelToSFF(outfile::Union{IOStream,IOBuffer}, model::AFeature, model_id::String,
                        targetloop::DNAString,
                        gene_exons::Dict{String,Int32},
                        maxlengths::Dict{String,Int32},
                        feature_stacks::AFeatureStack, strand::Char)
    gene = getFeatureName(first(model))
    expected_exons = gene_exons[gene]
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
        if gene ≠ "rps12B" && !isStartCodon(SubString(targetloop, first(model).start, first(model).start + 2), true, true)
            hasStart = false
        end
        stop_position = findfirst(isequal('*'), protein)
        if !isnothing(stop_position) && stop_position < length(protein)
            hasPrematureStop = true
        end
    end
    gene_length = last(model).start + last(model).length - first(model).start
    gene_length <= 0 && return
    for f in model
        stack = feature_stacks[findfirst(x->x.path == f.path, feature_stacks)]
        depth, coverage = getDepthAndCoverage(stack, f.start, f.length)
        depth == 0 && continue
        coverage == 0 && continue
        write(outfile, model_id)
        write(outfile, f.path[findall(isequal('/'), f.path)[2]:end])
        write(outfile, "\t")
        write(outfile, join([strand,string(f.start),string(f.length),string(f.phase)], "\t"))
        write(outfile, "\t")
        write(outfile, join([@sprintf("%.3f",f.length / stack.template.median_length),@sprintf("%.3f",depth),@sprintf("%.3f",coverage)], "\t"))
        write(outfile, "\t")
        if gene_length - get(maxlengths, gene, 0) < -10
            write(outfile, "possible pseudogene, shorter than 2nd copy ")
        end
        if exon_count < expected_exons
            write(outfile, "possible pseudogene, missing exon(s) ")
        end
        if !hasStart
            write(outfile, "possible pseudogene, no start codon ")
        end
        if hasPrematureStop
            write(outfile, "possible pseudogene, premature stop codon ")
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

function calc_maxlengths(fstrand_models::AAFeature,
    rstrand_models::AAFeature)::Dict{String,Int32}
    maxlengths = Dict{String,Int32}()
    for (fmodel, rmodel) in zip(fstrand_models, rstrand_models)
        if isempty(fmodel)
            @warn "forward model has no features"
        else
            gene = getFeatureName(first(fmodel))
            gene_length = last(fmodel).start + last(fmodel).length - first(fmodel).start
            maxlen = get(maxlengths, gene, 0)
            if gene_length > maxlen
                maxlengths[gene] = gene_length
            end
        end
        if isempty(rmodel)
            @warn "reverse model has no features"
        else
            gene = getFeatureName(first(rmodel))
            gene_length = last(rmodel).start + last(rmodel).length - first(rmodel).start
            maxlen = get(maxlengths, gene, 0)
            if gene_length > maxlen
                maxlengths[gene] = gene_length
            end
        end
    end
    maxlengths
end

function writeSFF(outfile::Union{String,IOStream,IOBuffer}, id::String, 
                fstrand_models::AAFeature,
                rstrand_models::AAFeature,
                gene_exons::Dict{String,Int32},
                fstrand_feature_stacks::AFeatureStack, 
                rstrand_feature_stacks::AFeatureStack, 
                targetloopf::DNAString, targetloopr::DNAString)

    maxlengths = calc_maxlengths(fstrand_models, rstrand_models)

    genome_length = length(fstrand_feature_stacks[1].stack)


    function out(outfile::Union{IOStream,IOBuffer})
        model_ids = Dict{String,Int32}()
        write(outfile, id, "\t", string(genome_length), "\n")
        for model in fstrand_models
            isempty(model) && continue
            model_id = getModelID!(model_ids, model)
            writeModelToSFF(outfile, model, model_id, targetloopf, gene_exons, maxlengths, fstrand_feature_stacks, '+')
        end
        for model in rstrand_models
            isempty(model) && continue
            model_id = getModelID!(model_ids, model)
            writeModelToSFF(outfile, model, model_id, targetloopr, gene_exons, maxlengths, rstrand_feature_stacks, '-')
        end
    end
    if typeof(outfile) == String
        maybe_gzopen(outfile::String, "w") do io
            out(io)
        end
    else
        out(outfile)
    end
end
