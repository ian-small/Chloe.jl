
include("Utilities.jl")
include("SuffixArrays.jl")
include("Alignments3.jl")
include("Annotations.jl")


using Dates

const ReferenceOrganisms = Dict(
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
    # 2* length(ReferenceOrganisms) from directory reference_1116
    refsrc::Array{String}
    refloops::Array{String}
    refSAs::Array{SuffixArray}
    refRAs::Array{SuffixArray}
    ref_features::Array{FeatureArray}
    # from .tsv file
    feature_templates::Array{FeatureTemplate}
    gene_exons::Dict{String,Int32}
end

function Base.show(io::IO, reference::Reference)
    t1 = "#templates=$(length(reference.feature_templates))"
    t2 = "#gene_exons=$(length(reference.gene_exons))[$(sum(values(reference.gene_exons)))]"
    t3 = "ref #seq=$(length(reference.refloops))[$(sum(map(x->length(x), reference.refloops)))]"
    print(io, "Reference: $t1, $t2, $t3")
end
"""
    readReferences(reference_dir, template_file_tsv)

creates a Reference object to feed into `annotate_one`
"""
function readReferences(refsdir::String, templates::String)::Reference

    num_refs = length(ReferenceOrganisms)

    refloops = Array{String}(undef, num_refs * 2)
    refSAs = Array{SuffixArray}(undef, num_refs * 2)
    refRAs = Array{SuffixArray}(undef, num_refs * 2)
    ref_features = Array{FeatureArray}(undef, num_refs * 2)
    refsrc = Array{String}(undef, num_refs * 2)
    if !isdir(refsdir)
        error("$(refsdir) is not a directory")
    end
    files = readdir(refsdir)
    idx = findall(x->endswith(x, ".sff"), files)
    if isempty(idx)
        error("no sff files found!")
    end
    iff_files = files[idx]

    for (i, ref) in enumerate(ReferenceOrganisms)
        path = joinpath(refsdir, ref.first * ".gwsas")
        if !isfile(path)
            error("no gwsas file $(path)")
        end
        refgwsas = readGenomeWithSAs(path, ref.first)
        rev = revComp(refgwsas.sequence)
        
        refloops[i * 2 - 1] = refgwsas.sequence * refgwsas.sequence[1:end - 1]
        refloops[i * 2] = rev * rev[1:end - 1]
        
        refSAs[i * 2 - 1] = refgwsas.forwardSA
        refSAs[i * 2] = refgwsas.reverseSA

        refRAs[i * 2 - 1] = makeSuffixArrayRanksArray(refgwsas.forwardSA)
        refRAs[i * 2] = makeSuffixArrayRanksArray(refgwsas.reverseSA)
        idx = findfirst(x->startswith(x, ref.second), iff_files)
        if idx === nothing
            error("no sff file for $(ref.second)")
        end
        feature_file = iff_files[idx]
        f_strand_features, r_strand_features = readFeatures(joinpath(refsdir, feature_file))
        
        ref_features[i * 2 - 1] = f_strand_features
        ref_features[i * 2] = r_strand_features

        refsrc[i * 2 - 1] = ref.first *  ":fwd"
        refsrc[i * 2] = ref.first * ":rev"

    end
    feature_templates, gene_exons = readTemplates(templates)
    return Reference(refsrc, refloops, refSAs, refRAs, ref_features, feature_templates, gene_exons)
end

const ns(td) = Time(Nanosecond(td))

MayBeString = Union{Nothing,String}
Strand = Tuple{AAFeature,AFeatureStack}

function do_strand(target_id::String, start_ns::UInt64, target_length::Int32,
    reference::Reference, coverages::Dict{String,Float32},
    strand::Char, blocks_aligned_to_target::Array{AlignedBlocks},
    targetloop::DNAString)::Strand

    annotations = Array{Annotation}(undef, 0)
    for (ref_feature_array, blocks) in zip(reference.ref_features, blocks_aligned_to_target)
        annotations = cat(annotations, findOverlaps(ref_feature_array, blocks), dims = 1)
    end
    sort!(annotations, by = x->x.path)
    strand_annotations = AnnotationArray(target_id, strand, annotations)

    @debug "[$(target_id)]$(strand) thread=$(Threads.threadid())"

    t4 = time_ns()
    @info "[$(target_id)]$(strand) overlapping ref annotations ($(length(annotations))): $(ns(t4 - start_ns))"

    strand_feature_stacks, shadow = fillFeatureStack(target_length, strand_annotations, reference.feature_templates)

    t5 = time_ns()
    @info "[$(target_id)]$(strand) ref features stacks ($(length(strand_feature_stacks))): $(ns(t5 - t4))"

    target_strand_features = FeatureArray(target_id, target_length, strand, AFeature(undef, 0))

    for stack in strand_feature_stacks
        left_border, length = alignTemplateToStack(stack, shadow)
        left_border == 0 && continue
        depth, coverage = getDepthAndCoverage(stack, left_border, length)
        if ((depth >= stack.template.threshold_counts) && (coverage >= stack.template.threshold_coverage))
            push!(target_strand_features.features, Feature(stack.path, left_border, length, 0))
        else
            @debug "[$(target_id)]$(strand) Below threshold: $(stack.path)"
        end
    end

    t6 = time_ns()
    @info "[$(target_id)]$(strand) aligning templates ($(length(target_strand_features.features))): $(ns(t6 - t5))"

    for feat in target_strand_features.features
        refineMatchBoundariesByOffsets!(feat, strand_annotations, target_length, coverages)
    end

    t7 = time_ns()
    @info "[$(target_id)]$(strand) refining match boundaries: $(ns(t7 - t6))"


    target_strand_models = groupFeaturesIntoGeneModels(target_strand_features)
    target_strand_models = refineGeneModels!(target_strand_models, target_length, targetloop, strand_annotations, strand_feature_stacks)

    t8 = time_ns()
    @info "[$(target_id)]$(strand) refining gene models: $(ns(t8 - t7))"
    return target_strand_models, strand_feature_stacks
end

"""
    annotate_one(references, fasta_file [,output_sff_file])

Annotate a single fasta file containting a *single* circular
DNA entry

writes an .sff file to `output_sff_file` or uses that id in the
fasta file to write `{target_id}.sff` in the current directory.

If output_sff_file is a *directory* write `{target_id}.sff` into that
directory.

`reference` are the reference annotations (see `readReferences`)
"""
MayBeIO = Union{String,IOBuffer,IOStream,Nothing}
function annotate_one(reference::Reference, fasta::Union{String,IOBuffer,IOStream},
    output::MayBeIO = nothing)

    num_refs = length(ReferenceOrganisms)
    t1 = time_ns()

    target_id, target_seqf = readFasta(fasta)
    target_length = Int32(length(target_seqf))
    
    @info "[$(target_id)] length: $(target_length)"
    
    target_seqr = revComp(target_seqf)
    targetloopf = target_seqf * target_seqf[1:end - 1]
    targetloopr = target_seqr * target_seqr[1:end - 1]
    
    target_saf = makeSuffixArray(targetloopf, true)
    target_raf = makeSuffixArrayRanksArray(target_saf)
    
    t2 = time_ns()

    @info "[$(target_id)] making suffix arrays: $(ns(t2 - t1)))"

    blocks_aligned_to_targetf = Array{AlignedBlocks}(undef, num_refs * 2)
    blocks_aligned_to_targetr = Array{AlignedBlocks}(undef, num_refs * 2)

    function alignit(refcount::Int)
        start = time_ns()
        refloop, refSA, refRA = reference.refloops[refcount], reference.refSAs[refcount], reference.refRAs[refcount]
        f_aligned_blocks, r_aligned_blocks = alignLoops(refloop, refSA, refRA, targetloopf, target_saf, target_raf)
        
        @debug "Coverage[$(Threads.threadid())][$(reference.refsrc[refcount])] ($(ns(time_ns() - start))): " forward = blockCoverage(f_aligned_blocks)  reverse = blockCoverage(r_aligned_blocks)

        # f_aligned_blocks contains matches between ref forward and target forward strands
        blocks_aligned_to_targetf[refcount] = f_aligned_blocks
        if refcount % 2 == 1 # aligning + strand to + strand
            blocks_aligned_to_targetr[refcount + 1] = r_aligned_blocks # r_aligned_blocks contains calculated matches between ref reverse and target reverse strands
        else    # aligning - strand to + strand
            blocks_aligned_to_targetr[refcount - 1] = r_aligned_blocks # r_aligned_blocks contains calculated matches between ref forward and target reverse strands
        end       
    end
    
    Threads.@threads for refno in 1:length(reference.refloops)
        alignit(refno)
    end

    t3 = time_ns()
    
    @info "[$(target_id)] aligning: ($(length(reference.refloops))) $(ns(t3 - t2))" 

    coverages = Dict{String,Float32}()
    for (i, ref) in enumerate(ReferenceOrganisms)
        coverage = 0
        coverage += blockCoverage(blocks_aligned_to_targetf[i * 2 - 1])
        coverage += blockCoverage(blocks_aligned_to_targetf[i * 2])
        coverages[ref[1]] = coverage /= target_length * 2
    end
    @debug "[$(target_id)] coverages:" coverages


    function watson(strands::Array{Strand})
        strands[1] = do_strand(target_id, t3, target_length, reference, coverages,
            '+', blocks_aligned_to_targetf, targetloopf)
    end
    function crick(strands::Array{Strand})
        strands[2] = do_strand(target_id, t3, target_length, reference, coverages,
            '-', blocks_aligned_to_targetr, targetloopr)
    end

    strands = Array{Strand}(undef, 2)
    Threads.@threads for worker in [watson, crick]
        worker(strands)
    end

    target_fstrand_models, fstrand_feature_stacks = strands[1]
    target_rstrand_models, rstrand_feature_stacks = strands[2]

    if output != nothing
        if typeof(output) == String
            fname = output::String
            if isdir(fname)
                fname = joinpath(fname, "$(target_id).sff")
            end
        else
            fname = output
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
function annotate_one(reference::Reference, fasta::Union{String,IOBuffer,IOStream})
    annotate_one(reference, fasta, IOBuffer())
end
function annotate(refsdir::String, templates::String, fa_files::Array{String}, output::MayBeString)

    reference = readReferences(refsdir, templates)

    for infile in fa_files
        annotate_one(reference, infile, output)
    end
end
