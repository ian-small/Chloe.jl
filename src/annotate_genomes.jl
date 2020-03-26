

include("Utilities.jl")
include("SuffixArrays.jl")
include("Annotations.jl")
include("Alignments3.jl")


using Dates

# refsdir = ARGS[1]
# const refs = Dict("NC_004543"=>"Anthoceros","AP000423"=>"Arabidopsis","MF177093"=>"Azolla","NC_001319"=>"Marchantia","NC_022137"=>"Marsilea",
# "JX512022"=>"Medicago","Z00044"=>"Nicotiana","AP005672"=>"Physcomitrella","KT634228"=>"Picea","FJ755183"=>"Selaginella","NC_002202"=>"Spinacia","NC_001666"=>"Zea")
# const refs = Dict("AP000423"=>"Arabidopsis","JX512022"=>"Medicago","Z00044"=>"Nicotiana","KT634228"=>"Picea","NC_002202"=>"Spinacia","NC_001666"=>"Zea")
# const refs = Dict("NC_001666"=>"Zea")
# const refs = Dict("BK010421"=>"Arabidopsis")
const refs = Dict(
    "AP000423"  => "Arabidopsis",
    "JX512022"  => "Medicago",
    "Z00044"    => "Nicotiana",
    "KT634228"  => "Picea",
    "NC_002202" => "Spinacia",
    "NC_001666" => "Zea",
    "NC_005086" => "Amborella",
    "NC_016986" => "Ginkgo",
    "NC_021438" => "Gnetum",
    "NC_030504" => "Liriodendron",
    "NC_024542" => "Nymphaea",
    "NC_031333" => "Oryza",
    "NC_026040" => "Zamia"
)

struct Reference
    refloops::Vector{String}
    refSAs::Array{Array{Int32}}
    refRAs::Array{Array{Int32}}
    ref_features::Array{FeatureArray}
    feature_templates::Array{FeatureTemplate}
    gene_exons::Dict{String,Int32}
end

function show_reference(reference::Reference)
    return """templates=$(length(reference.feature_templates)), 
        gene_exons=$(length(reference.gene_exons))[$(sum(values(reference.gene_exons)))],
        ref seq loops=$(length(reference.refloops))[$(sum(map(x->length(x), reference.refloops)))]
        """
end

function readReferences(refsdir::String, templates::String)

    num_refs = length(refs)

    refloops = Vector{String}(undef, num_refs * 2)
    refSAs = Array{Array{Int32,1}}(undef, num_refs * 2)
    refRAs = Array{Array{Int32,1}}(undef, num_refs * 2)
    ref_features = Array{FeatureArray,1}(undef, num_refs * 2)

    files = readdir(refsdir)
    iff_files = files[findall(x->endswith(x, ".sff"), files)]

    for (i, ref) in enumerate(refs)
        refgwsas = readGenomeWithSAs(joinpath(refsdir, ref.first * ".gwsas"), ref.first)
        refloops[i * 2 - 1] = refgwsas.sequence * refgwsas.sequence[1:end - 1]
        rev = revComp(refgwsas.sequence)
        refloops[i * 2] = rev * rev[1:end - 1]
        refSAs[i * 2 - 1] = refgwsas.forwardSA
        refSAs[i * 2] = refgwsas.reverseSA

        refRAs[i * 2 - 1] = makeSuffixArrayRanksArray(refgwsas.forwardSA)
        refRAs[i * 2] = makeSuffixArrayRanksArray(refgwsas.reverseSA)

        # feature_file = iff_files[findfirst(x->startswith(x,ref.second),iff_files)]
        feature_file = iff_files[findfirst(x->startswith(x, ref.second), iff_files)]
        f_strand_features, r_strand_features = readFeatures(joinpath(refsdir, feature_file))
        ref_features[i * 2 - 1] = f_strand_features
        ref_features[i * 2] = r_strand_features
    end
    feature_templates, gene_exons = readTemplates(templates)
    return Reference(refloops, refSAs, refRAs, ref_features, feature_templates, gene_exons)
end

const ns(td) = Time(Nanosecond(td))

MayBeString = Union{Nothing,String}

function annotate_one(fasta::String, reference::Reference, output::MayBeString)

    num_refs = length(refs)
    t1 = time_ns()
    if !isfile(fasta)
        error("$(fasta): not a file!")
    end
    target_id, target_seqf = readFasta(fasta)
    target_length = length(target_seqf)
    
    @info "[$(target_id)] length: $(target_length)"
    
    target_seqr = revComp(target_seqf)
    targetloopf = target_seqf * target_seqf[1:end - 1]
    target_saf = makeSuffixArray(targetloopf, true)
    target_raf = makeSuffixArrayRanksArray(target_saf)

    t2 = time_ns()

    @info "[$(target_id)] making suffix arrays: $(ns(t2 - t1)))"

    blocks_aligned_to_targetf = Array{Array{Tuple{Int32,Int32,Int32}}}(undef, num_refs * 2)
    blocks_aligned_to_targetr = Array{Array{Tuple{Int32,Int32,Int32}}}(undef, num_refs * 2)
    # refcount = 0
    # for (refcount, (refloop, refSA, refRA)) in enumerate(zip(reference.refloops, reference.refSAs, reference.refRAs))
    #     # refcount += 1
    #     f_aligned_blocks, r_aligned_blocks = alignLoops(refloop, refSA, refRA, targetloopf, target_saf, target_raf)
    #     if refcount % 2 == 1 # aligning + strand to + strand
    #         blocks_aligned_to_targetf[refcount] = f_aligned_blocks # f_aligned_blocks contains matches between ref forward and target forward strands
    #         blocks_aligned_to_targetr[refcount + 1] = r_aligned_blocks # r_aligned_blocks contains calculated matches between ref reverse and target reverse strands
    #     else    # aligning - strand to + strand
    #         blocks_aligned_to_targetf[refcount] = f_aligned_blocks # f_aligned_blocks contains matches between ref reverse and target forward strands
    #         blocks_aligned_to_targetr[refcount - 1] = r_aligned_blocks # r_aligned_blocks contains calculated matches between ref forward and target reverse strands
    #     end
    # end
    function alignit(refcount)
        refloop, refSA, refRA = reference.refloops[refcount], reference.refSAs[refcount], reference.refRAs[refcount]
        f_aligned_blocks, r_aligned_blocks = alignLoops(refloop, refSA, refRA, targetloopf, target_saf, target_raf)
        if refcount % 2 == 1 # aligning + strand to + strand
            blocks_aligned_to_targetf[refcount] = f_aligned_blocks # f_aligned_blocks contains matches between ref forward and target forward strands
            blocks_aligned_to_targetr[refcount + 1] = r_aligned_blocks # r_aligned_blocks contains calculated matches between ref reverse and target reverse strands
        else    # aligning - strand to + strand
            blocks_aligned_to_targetf[refcount] = f_aligned_blocks # f_aligned_blocks contains matches between ref reverse and target forward strands
            blocks_aligned_to_targetr[refcount - 1] = r_aligned_blocks # r_aligned_blocks contains calculated matches between ref forward and target reverse strands
        end       
    end
    Threads.@threads for refno in 1:length(reference.refloops) 
        alignit(refno)
    end

    t3 = time_ns()
    
    @info "[$(target_id)] aligning: $(ns(t3 - t2))" 

    coverages = Dict{String,Real}()
    for (i, ref) in enumerate(refs)
        coverage = 0
        coverage += blockCoverage(blocks_aligned_to_targetf[i * 2 - 1])
        coverage += blockCoverage(blocks_aligned_to_targetf[i * 2])
        coverages[ref[1]] = coverage /= target_length * 2
    end
    @debug "[$(target_id)] coverages:" coverages


    # f_strand_annotations = AnnotationArray(target_id, '+', Array{Annotation,1}(undef, 0))
    # r_strand_annotations = AnnotationArray(target_id, '-', Array{Annotation,1}(undef, 0))
    # for (ref_feature_array, blocksf, blocksr) in zip(reference.ref_features, blocks_aligned_to_targetf, blocks_aligned_to_targetr)
    #     f_strand_annotations.annotations = cat(f_strand_annotations.annotations, pushFeatures(ref_feature_array, target_id, '+', blocksf).annotations, dims = 1)
    #     r_strand_annotations.annotations = cat(r_strand_annotations.annotations, pushFeatures(ref_feature_array, target_id, '-', blocksr).annotations, dims = 1)
    # end
    # sort!(f_strand_annotations.annotations, by = x->x.path)
    # sort!(r_strand_annotations.annotations, by = x->x.path)

    f_annotations = Array{Annotation,1}(undef, 0)
    r_annotations = Array{Annotation,1}(undef, 0)
    for (ref_feature_array, blocksf, blocksr) in zip(reference.ref_features, blocks_aligned_to_targetf, blocks_aligned_to_targetr)
        f_annotations = cat(f_annotations, pushFeatures(ref_feature_array, target_id, '+', blocksf).annotations, dims = 1)
        r_annotations = cat(r_annotations, pushFeatures(ref_feature_array, target_id, '-', blocksr).annotations, dims = 1)
    end
    sort!(f_annotations, by = x->x.path)
    sort!(r_annotations, by = x->x.path)
    f_strand_annotations = AnnotationArray(target_id, '+', f_annotations)
    r_strand_annotations = AnnotationArray(target_id, '-', r_annotations)

    t4 = time_ns()
    @info "[$(target_id)] pushing annotations: $(ns(t4 - t3))"

    fstrand_feature_stacks, fshadow = stackFeatures(length(target_seqf), f_strand_annotations, reference.feature_templates)
    rstrand_feature_stacks, rshadow = stackFeatures(length(target_seqr), r_strand_annotations, reference.feature_templates)

    t5 = time_ns()
    @info "[$(target_id)] stacking features: $(ns(t5 - t3))"

    target_fstrand_features = FeatureArray(target_id, target_length, '+', Array{Feature,1}(undef, 0))

    for (i, fstack) in enumerate(fstrand_feature_stacks)
        left_border, length = alignTemplateToStack(fstack, fshadow)
        left_border == 0 && continue
        depth, coverage = getDepthAndCoverage(fstack, left_border, length)
            # println(fstack.path," ",depth," ",coverage)
        if ((depth >= fstack.template.threshold_counts) && (coverage >= fstack.template.threshold_coverage))
            push!(target_fstrand_features.features, Feature(fstack.path, left_border, length, 0))
        else
            # println("Below threshold: " * fstack.path)
        end
    end

    target_rstrand_features = FeatureArray(target_id, target_length, '-', Array{Feature,1}(undef, 0))

    for (i, rstack) in enumerate(rstrand_feature_stacks)
        left_border, length = alignTemplateToStack(rstack, rshadow)
        left_border == 0 && continue
        depth, coverage = getDepthAndCoverage(rstack, left_border, length)
        if ((depth >= rstack.template.threshold_counts) && (coverage >= rstack.template.threshold_coverage))
            push!(target_rstrand_features.features, Feature(rstack.path, left_border, length, 0))
        else
            # println("Below threshold: " * fstack.path)
        end
    end

    t6 = time_ns()
    @info "[$(target_id)] aligning templates: $(ns(t6 - t5))"

    for feat in target_fstrand_features.features
        feat = refineMatchBoundariesByOffsets!(feat, f_strand_annotations, target_length, coverages)
    end

    for feat in target_rstrand_features.features
        feat = refineMatchBoundariesByOffsets!(feat, r_strand_annotations, target_length, coverages)
    end

    t7 = time_ns()
    @info "[$(target_id)] refining match boundaries: $(ns(t7 - t6))"

    target_fstrand_models = groupFeaturesIntoGeneModels(target_fstrand_features)
    target_fstrand_models = refineGeneModels!(target_length, targetloopf, target_fstrand_models, f_strand_annotations, fstrand_feature_stacks)
    target_rstrand_models = groupFeaturesIntoGeneModels(target_rstrand_features)
    targetloopr = target_seqr * target_seqr[1:end - 1]
    target_rstrand_models = refineGeneModels!(target_length, targetloopr, target_rstrand_models, r_strand_annotations, rstrand_feature_stacks)

    t8 = time_ns()
    @info "[$(target_id)] refining gene models: $(ns(t8 - t7))"

    if output != nothing
        fname = output::String
        if isdir(fname)
            fname = joinpath(fname, "$(target_id).sff")
        end
    else
        fname = "$(target_id).sff"
    end
    writeSFF(fname, target_id, target_fstrand_models, target_rstrand_models,
        reference.gene_exons, fstrand_feature_stacks, rstrand_feature_stacks,
        targetloopf, targetloopr)
    @info "[$(target_id)] Overall: $(ns(time_ns() - t1))"
    return fname, target_id

end

function annotate(refsdir::String, templates::String, fa_files::Array{String,1}, output::MayBeString)

    reference = readReferences(refsdir, templates)

    for infile in fa_files
        annotate_one(infile, reference, output)
    end
end
