include("Utilities.jl")
include("SuffixArrays.jl")
include("Annotations.jl")
include("Alignments3.jl")

refsdir = ARGS[1]
#const refs = Dict("NC_004543"=>"Anthoceros","AP000423"=>"Arabidopsis","MF177093"=>"Azolla","NC_001319"=>"Marchantia","NC_022137"=>"Marsilea","JX512022"=>"Medicago","Z00044"=>"Nicotiana","AP005672"=>"Physcomitrella","KT634228"=>"Picea","FJ755183"=>"Selaginella","NC_002202"=>"Spinacia","NC_001666"=>"Zea")
#const refs = Dict("AP000423"=>"Arabidopsis","JX512022"=>"Medicago","Z00044"=>"Nicotiana","KT634228"=>"Picea","NC_002202"=>"Spinacia","NC_001666"=>"Zea")
#const refs = Dict("NC_001666"=>"Zea")
#const refs = Dict("BK010421"=>"Arabidopsis")
const refs = Dict("AP000423"=>"Arabidopsis","JX512022"=>"Medicago","Z00044"=>"Nicotiana","KT634228"=>"Picea",
    "NC_002202"=>"Spinacia","NC_001666"=>"Zea","NC_005086"=>"Amborella","NC_016986"=>"Ginkgo","NC_021438"=>"Gnetum",
    "NC_030504"=>"Liriodendron","NC_024542"=>"Nymphaea","NC_031333"=>"Oryza","NC_026040"=>"Zamia")

num_refs = length(refs)

refloops = Vector{String}(undef,num_refs*2)
refSAs = Array{Array{Int32,1}}(undef,num_refs*2)
refRAs = Array{Array{Int32,1}}(undef,num_refs*2)
ref_features = Array{FeatureArray,1}(undef,num_refs*2)

files = readdir(refsdir)
iff_files = files[findall(x->endswith(x,".sff"), files)]

for (i,ref) in enumerate(refs)
    refgwsas = readGenomeWithSAs(joinpath(refsdir,ref.first*".gwsas"),ref.first)
    refloops[i*2-1] = refgwsas.sequence * refgwsas.sequence[1:end-1]
    rev = revComp(refgwsas.sequence)
    refloops[i*2] = rev * rev[1:end-1]
    refSAs[i*2-1] = refgwsas.forwardSA
    refSAs[i*2] = refgwsas.reverseSA

    refRAs[i*2-1] = makeSuffixArrayRanksArray(refgwsas.forwardSA)
    refRAs[i*2] = makeSuffixArrayRanksArray(refgwsas.reverseSA)

    #feature_file = iff_files[findfirst(x->startswith(x,ref.second),iff_files)]
    feature_file = iff_files[findfirst(x->startswith(x,ref.second),iff_files)]
    f_strand_features,r_strand_features = readFeatures(joinpath(refsdir,feature_file))
    ref_features[i*2-1] = f_strand_features
    ref_features[i*2] = r_strand_features
end

const feature_templates,gene_exons = readTemplates(ARGS[2])

for infile in ARGS[3:end]

    t1 = time_ns()
    target_id,target_seqf = readFasta(infile)
    target_length = length(target_seqf)
    println(target_id)
    target_seqr = revComp(target_seqf)
    targetloopf = target_seqf*target_seqf[1:end-1]
    target_saf = makeSuffixArray(targetloopf,true)
    target_raf = makeSuffixArrayRanksArray(target_saf)

    t2 = time_ns()
    #println("making suffix arrays: " * string(t2-t1))

    blocks_aligned_to_targetf = Array{Array{Tuple{Int32, Int32, Int32}}}(undef,num_refs*2)
    blocks_aligned_to_targetr = Array{Array{Tuple{Int32, Int32, Int32}}}(undef,num_refs*2)
    refcount = 0
    for (refloop,refSA,refRA) in zip(refloops,refSAs,refRAs)
        refcount += 1
        f_aligned_blocks,r_aligned_blocks = alignLoops(refloop,refSA,refRA,targetloopf,target_saf,target_raf)
        if refcount % 2 == 1 # aligning + strand to + strand
            blocks_aligned_to_targetf[refcount] = f_aligned_blocks # f_aligned_blocks contains matches between ref forward and target forward strands
            blocks_aligned_to_targetr[refcount+1] = r_aligned_blocks # r_aligned_blocks contains calculated matches between ref reverse and target reverse strands
        else    # aligning - strand to + strand
            blocks_aligned_to_targetf[refcount] = f_aligned_blocks # f_aligned_blocks contains matches between ref reverse and target forward strands
            blocks_aligned_to_targetr[refcount-1] = r_aligned_blocks # r_aligned_blocks contains calculated matches between ref forward and target reverse strands
        end
    end

    coverages = Dict{String,Real}()
    for (i,ref) in enumerate(refs)
        coverage = 0
        coverage += blockCoverage(blocks_aligned_to_targetf[i*2-1])
        coverage += blockCoverage(blocks_aligned_to_targetf[i*2])
        coverages[ref[1]] = coverage /= target_length*2
    end
    #println(coverages)

    t3 = time_ns()
    #println("aligning: " * string(t3-t2))

    f_strand_annotations = AnnotationArray(target_id,'+',Array{Annotation,1}(undef,0))
    r_strand_annotations = AnnotationArray(target_id,'-',Array{Annotation,1}(undef,0))
    for (ref_feature_array,blocksf,blocksr) in zip(ref_features,blocks_aligned_to_targetf,blocks_aligned_to_targetr)
        f_strand_annotations.annotations = cat(f_strand_annotations.annotations,pushFeatures(ref_feature_array, target_id, '+', blocksf).annotations,dims=1)
        r_strand_annotations.annotations = cat(r_strand_annotations.annotations,pushFeatures(ref_feature_array, target_id, '-', blocksr).annotations,dims=1)
    end
    sort!(f_strand_annotations.annotations, by = x -> x.path)
    sort!(r_strand_annotations.annotations, by = x -> x.path)

    t4 = time_ns()
    #println("pushing annotations: " * string(t4-t3))

    fstrand_feature_stacks, fshadow = stackFeatures(length(target_seqf),f_strand_annotations,feature_templates)
    rstrand_feature_stacks, rshadow = stackFeatures(length(target_seqr),r_strand_annotations,feature_templates)

    t5 = time_ns()
    #println("stacking features: " * string(t5-t4))

    target_fstrand_features=FeatureArray(target_id,target_length,'+',Array{Feature,1}(undef,0))

    for (i,fstack) in enumerate(fstrand_feature_stacks)
        left_border, length = alignTemplateToStack(fstack, fshadow)
        left_border == 0 && continue
        depth, coverage = getDepthAndCoverage(fstack, left_border, length)
            #println(fstack.path," ",depth," ",coverage)
        if ((depth >= fstack.template.threshold_counts) && (coverage >= fstack.template.threshold_coverage))
            push!(target_fstrand_features.features,Feature(fstack.path, left_border, length,0))
        else
            #println("Below threshold: " * fstack.path)
        end
    end

    target_rstrand_features=FeatureArray(target_id,target_length,'-',Array{Feature,1}(undef,0))

    for (i,rstack) in enumerate(rstrand_feature_stacks)
        left_border, length = alignTemplateToStack(rstack, rshadow)
        left_border == 0 && continue
        depth, coverage = getDepthAndCoverage(rstack, left_border, length)
        if ((depth >= rstack.template.threshold_counts) && (coverage >= rstack.template.threshold_coverage))
            push!(target_rstrand_features.features,Feature(rstack.path, left_border, length,0))
        else
            #println("Below threshold: " * fstack.path)
        end
    end

    t6 = time_ns()
    #println("aligning templates: " * string(t6-t5))

    for feat in target_fstrand_features.features
        feat = refineMatchBoundariesByOffsets!(feat,f_strand_annotations,target_length,coverages)
    end

    for feat in target_rstrand_features.features
        feat = refineMatchBoundariesByOffsets!(feat,r_strand_annotations,target_length,coverages)
    end

    t7 = time_ns()
    #println("refining match boundaries: " * string(t7-t6))

    target_fstrand_models = groupFeaturesIntoGeneModels(target_fstrand_features)
    target_fstrand_models = refineGeneModels!(target_length,targetloopf,target_fstrand_models,f_strand_annotations,fstrand_feature_stacks)
    target_rstrand_models = groupFeaturesIntoGeneModels(target_rstrand_features)
    targetloopr = target_seqr*target_seqr[1:end-1]
    target_rstrand_models = refineGeneModels!(target_length,targetloopr,target_rstrand_models,r_strand_annotations,rstrand_feature_stacks)

    t8 = time_ns()
    #println("refining gene models: " * string(t8-t7))
    #println("Overall: " * string(t8-t1))

    writeSFF(target_id*".sff",target_id,target_fstrand_models,target_rstrand_models,gene_exons,fstrand_feature_stacks,rstrand_feature_stacks,targetloopf,targetloopr)
end
