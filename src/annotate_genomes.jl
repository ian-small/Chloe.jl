DNAString = AbstractString

include("MMappedString.jl")
include("Utilities.jl")
include("SuffixArrays.jl")
include("Alignments3.jl")
include("Annotations.jl")


# can be any of String, MMappedString, MappedPtrString{ASCII}
# If you have memory mapped files (see julia chloe.jl mmap *.fa)
# then use MMappedString{ASCII} since it will not read
# the data backing the memory mapping
MappedPtrString = MMappedString{ASCII}
# MappedPtrString = String


# import Dates: Time, Nanosecond
import Base
import Printf: @sprintf
import JSON
import Mmap


# const ReferenceOrganisms = Dict(
#     "AP000423"  => "Arabidopsis",
#     "JX512022"  => "Medicago",
#     "Z00044"    => "Nicotiana",
#     "KT634228"  => "Picea",
#     "NC_002202" => "Spinacia",
#     "NC_001666" => "Zea",
#     "NC_005086" => "Amborella",
#     "NC_016986" => "Ginkgo",
#     "NC_021438" => "Gnetum",
#     "NC_030504" => "Liriodendron",
#     "NC_024542" => "Nymphaea",
#     "NC_031333" => "Oryza",
#     "NC_026040" => "Zamia"
# )

struct FwdRev{T}
    forward::T
    reverse::T
end

Base.:(==)(x::FwdRev{T}, y::FwdRev{T}) where T = x.forward == y.forward && x.reverse == y.reverse

function mmap_suffix_arrays(filename::String)
    size = filesize(filename)
    nelem = floor(Int, (size + 2) / 20)
    
    nbytes = 4 * nelem
    sbytes = 2 * nelem - 1

    f = open(filename)
    # zero offset???
    saf = Mmap.mmap(f, SuffixArray, nelem, 0 + 0 * nbytes)
    sar = Mmap.mmap(f, SuffixArray, nelem, 0 + 1 * nbytes)
    raf = Mmap.mmap(f, SuffixArray, nelem, 0 + 2 * nbytes)
    rar = Mmap.mmap(f, SuffixArray, nelem, 0 + 3 * nbytes)
    sf =  Mmap.mmap(f, Vector{UInt8}, sbytes, 0 + 4 * nbytes)
    sr  = Mmap.mmap(f, Vector{UInt8}, sbytes, 0 + 4 * nbytes + sbytes)
    msf, msr = MappedPtrString(sf),  MappedPtrString(sr)
    return FwdRev(saf, sar), FwdRev(raf, rar), FwdRev(msf, msr)
end
struct Reference
    # length(ReferenceOrganisms) from directory reference_1116
    refsrc::Vector{String}
    refloops::Vector{FwdRev{MappedPtrString}}
    refSAs::Vector{FwdRev{SuffixArray}}
    refRAs::Vector{FwdRev{SuffixArray}}
    ref_features::Vector{FwdRev{FeatureArray}}
    # any open memory mapped files...
    # _mmaps::Vector{Union{Nothing,IOStream}}
    # from .tsv file
    feature_templates::Vector{FeatureTemplate}
    gene_exons::Dict{String,Int32}
    # referenceOrganisms::Dict{String,String}
end

datasize(t::T) where T = sizeof(t)
datasize(f::FwdRev{T}) where T = sizeof(FwdRev{T}) + datasize(f.forward) + datasize(f.reverse)
datasize(v::Vector{T}) where T = sum(datasize(a) for a in v)
datasize(t::Dict{K,V}) where {K,V} = sum(datasize(e.first) + datasize(e.second) for e in t)
datasize(m::MMappedString) = length(m.ptr)
datasize(r::Reference) = begin
    (sizeof(Reference)
    + datasize(r.refsrc) 
    + datasize(r.refloops)
    + datasize(r.refSAs)
    + datasize(r.refRAs)
    + datasize(r.ref_features)
    + datasize(r.feature_templates)
    + datasize(r.gene_exons))
end


function human(num::Int)::String
    if num === 0
        return "0B"
    end
    magnitude = floor(Int, log10(abs(num)) / 3)
    val = num / (1000^magnitude)
    sval = @sprintf("%.1f", val)
    if magnitude > 7
        return "$(sval)YB"
    end
    p = ["", "k", "M", "G", "T", "P", "E", "Z"][magnitude + 1]
    return "$(sval)$(p)B"
end
# stops the REPL printing the entire 7MB of sequences!
function Base.show(io::IO, reference::Reference)
    function fr_length(fr)
        return length(fr.forward) + length(fr.reverse)
    end

    sabytes = 2 * sizeof(Int32) * (reference.refSAs .|> fr_length |> sum)
    bp = reference.refloops .|> fr_length |> sum
    
    t1 = "#templates=$(reference.feature_templates |> length)"
    t2 = "#gene_exons=$(reference.gene_exons |> length)[$(reference.gene_exons |> values |> sum)]"
    t3 = "#seq=$(2 * (reference.refloops |> length))[$(human(bp)) bp]"
    print(io, "Reference: $t1, $t2, $t3, suffix=$(human(sabytes)), total=$(human(sabytes + bp))")

end
   


function read_gwsas(gwsas::String, id::String)
    
    function FwdRevSA(gwsas::GenomeWithSAs)
        FwdRev(gwsas.forwardSA, gwsas.reverseSA)
    end
    function FwdRevRA(gwsas::GenomeWithSAs)
        FwdRev(makeSuffixArrayRanksArray(gwsas.forwardSA), makeSuffixArrayRanksArray(gwsas.reverseSA))
    end
    function FwdRevMap(fwd::String, rev::String)
        # makes a copy unfortuately...
        FwdRev(MappedPtrString(fwd), MappedPtrString(rev))
    end
         
    refgwsas = readGenomeWithSAs(gwsas, id)
    fwd = refgwsas.sequence
    rev = revComp(fwd)
    fwd = fwd * fwd[1:end - 1]
    rev = rev * rev[1:end - 1]


    FwdRevSA(refgwsas), FwdRevRA(refgwsas), FwdRevMap(fwd, rev)
end
"""
    readReferences(reference_dir, template_file_tsv)

creates a Reference object to feed into `annotate_one`
"""
function readReferences(refsdir::String, templates::String)::Reference

    if !isdir(refsdir)
        error("$(refsdir) is not a directory")
    end
    ReferenceOrganisms = open(joinpath(refsdir, "ReferenceOrganisms.json")) do f
        JSON.parse(f, dicttype=Dict{String,String})
    end

    num_refs = length(ReferenceOrganisms)

    refloops = Vector{FwdRev{MappedPtrString}}(undef, num_refs)
    refSAs = Vector{FwdRev{SuffixArray}}(undef, num_refs)
    refRAs = Vector{FwdRev{SuffixArray}}(undef, num_refs)
    ref_features = Vector{FwdRev{FeatureArray}}(undef, num_refs)
    refsrc = Vector{String}(undef, num_refs)
    # mmapped = Vector{Union{Nothing,IOStream}}(undef, num_refs)
    

    files = readdir(refsdir)
    idx = findall(x -> endswith(x, ".sff"), files)
    if isempty(idx)
        error("No sff files found!")
    end
    
    sff_files = files[idx]

    for (i, ref) in enumerate(ReferenceOrganisms)
        gwsas = joinpath(refsdir, ref.first * ".gwsas")
        mmap = joinpath(refsdir, ref.first * ".mmap")
        
        if isfile(mmap)
            refSAs[i], refRAs[i], refloops[i] = mmap_suffix_arrays(mmap)
            @info "found mmap file for: $(ref.first)"
        elseif isfile(gwsas)
            refSAs[i], refRAs[i], refloops[i] = read_gwsas(gwsas, ref.first)
            @info "found gwsas file for: $(ref.first)"
        else
            error("no data file for $(ref.first)")
        end
        idx = findfirst(x -> startswith(x, ref.second), sff_files)
        if idx === nothing
            # TODO check for fasta files
            error("no sff file for $(ref.first) -> $(ref.second)")
        end
        feature_file = sff_files[idx]
        f_strand_features, r_strand_features = readFeatures(joinpath(refsdir, feature_file))

        ref_features[i] = FwdRev(f_strand_features, r_strand_features)
        refsrc[i] = ref.first


    end
    feature_templates, gene_exons = readTemplates(templates)
    return Reference(refsrc, refloops, refSAs, refRAs, ref_features, feature_templates, gene_exons)
end




# const ns(td) = Time(Nanosecond(td))
const ns(td) = @sprintf("%.3fs", td / 1e9)

MayBeString = Union{Nothing,String}
Strand = Tuple{AAFeature,AFeatureStack}
AAlignedBlocks = Vector{FwdRev{AlignedBlocks}}

function do_strand(target_id::String, start_ns::UInt64, target_length::Int32,
    reference::Reference, coverages::Dict{String,Float32},
    strand::Char, blocks_aligned_to_target::AAlignedBlocks,
    targetloop::DNAString)::Strand

    annotations = Vector{Annotation}()
    for (ref_feature_array, blocks) in zip(reference.ref_features, blocks_aligned_to_target)
        annotations = cat(annotations, findOverlaps(ref_feature_array.forward, blocks.forward), dims=1)
        annotations = cat(annotations, findOverlaps(ref_feature_array.reverse, blocks.reverse), dims=1)
    end
    sort!(annotations, by=x -> x.path)
    strand_annotations = AnnotationArray(target_id, strand, annotations)

    @debug "[$target_id]$strand thread=$(Threads.threadid())"

    t4 = time_ns()
    @info "[$target_id]$strand overlapping ref annotations ($(length(annotations))): $(ns(t4 - start_ns))"

    strand_feature_stacks, shadow = fillFeatureStack(target_length, strand_annotations, reference.feature_templates)

    t5 = time_ns()
    @info "[$target_id]$strand ref features stacks ($(length(strand_feature_stacks))): $(ns(t5 - t4))"

    target_strand_features = FeatureArray(target_id, target_length, strand, AFeature())

    for stack in strand_feature_stacks
        left_border, length = alignTemplateToStack(stack, shadow)
        left_border == 0 && continue
        depth, coverage = getDepthAndCoverage(stack, left_border, length)
        if ((depth >= stack.template.threshold_counts) && (coverage >= stack.template.threshold_coverage))
            push!(target_strand_features.features, Feature(stack.path, left_border, length, 0))
        else
            @debug "[$target_id]$strand Below threshold: $(stack.path)"
        end
    end

    t6 = time_ns()
    @info "[$target_id]$strand aligning templates ($(length(target_strand_features.features))): $(ns(t6 - t5))"

    for feat in target_strand_features.features
        refineMatchBoundariesByOffsets!(feat, strand_annotations, target_length, coverages)
    end

    t7 = time_ns()
    @info "[$target_id]$strand refining match boundaries: $(ns(t7 - t6))"


    target_strand_models = groupFeaturesIntoGeneModels(target_strand_features)
    target_strand_models = refineGeneModels!(target_strand_models, target_length, targetloop, strand_annotations, strand_feature_stacks)

    t8 = time_ns()
    @info "[$target_id]$strand refining gene models: $(ns(t8 - t7))"
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
MayBeIO = Union{String,IO,Nothing}
function annotate_one(reference::Reference, fasta::Union{String,IO},
    output::MayBeIO=nothing)

    num_refs = length(reference.refsrc)
    t1 = time_ns()

    target_id, target_seqf = readFasta(fasta)
    target_length = Int32(length(target_seqf))
    
    @info "[$target_id] length: $target_length"
    
    target_seqr = revComp(target_seqf)
    
    targetloopf = target_seqf * target_seqf[1:end - 1]
    target_saf = makeSuffixArray(targetloopf, true)
    target_raf = makeSuffixArrayRanksArray(target_saf)

    targetloopr = target_seqr * target_seqr[1:end - 1]
    target_sar = makeSuffixArray(targetloopr, true)
    target_rar = makeSuffixArrayRanksArray(target_sar)
    
    targetloopf = MappedPtrString(targetloopf)
    targetloopr = MappedPtrString(targetloopr)

    t2 = time_ns()

    @info "[$target_id] made suffix arrays: $(ns(t2 - t1))"

    blocks_aligned_to_targetf = AAlignedBlocks(undef, num_refs)
    blocks_aligned_to_targetr = AAlignedBlocks(undef, num_refs)

    function alignit(refcount::Int)
        start = time_ns()
        refloop, refSA, refRA = reference.refloops[refcount], reference.refSAs[refcount], reference.refRAs[refcount]
        
        ff, fr = alignLoops(refloop.forward, refSA.forward, refRA.forward, targetloopf, target_saf, target_raf) 
        rf, rr = alignLoops(refloop.reverse, refSA.reverse, refRA.reverse, targetloopf, target_saf, target_raf)
        # *OR* use only forward
        # rr_aligned_blocks, rf_aligned_blocks = alignLoops(refloop.forward, refSA.forward, refRA.forward, targetloopr, target_sar, target_rar)
        # note cross ...
        blocks_aligned_to_targetf[refcount] = FwdRev(ff, rf)
        blocks_aligned_to_targetr[refcount] = FwdRev(rr, fr)
        
        @info "[$target_id]± aligned $(reference.refsrc[refcount]) ($(length(ff)),$(length(rf))) $(ns(time_ns() - start))"

    end

    Threads.@threads for refno in 1:length(reference.refloops)
        alignit(refno)
    end

    t3 = time_ns()
    
    @info "[$target_id] aligned: ($(length(reference.refloops))) $(ns(t3 - t2))" 

    coverages = Dict{String,Float32}()
    for (i, ref) in enumerate(reference.refsrc)
        coverage = 0
        coverage += blockCoverage(blocks_aligned_to_targetf[i].forward)
        coverage += blockCoverage(blocks_aligned_to_targetf[i].reverse)
        coverages[ref] = coverage /= target_length * 2
    end
    @debug "[$target_id] coverages:" coverages


    function watson(strands::Vector{Strand})
        strands[1] = do_strand(target_id, t3, target_length, reference, coverages,
            '+', blocks_aligned_to_targetf, targetloopf)
    end
    function crick(strands::Vector{Strand})
        strands[2] = do_strand(target_id, t3, target_length, reference, coverages,
            '-', blocks_aligned_to_targetr, targetloopr)
    end

    strands = Vector{Strand}(undef, 2)
    Threads.@threads for worker in [watson, crick]
        worker(strands)
    end

    target_fstrand_models, fstrand_feature_stacks = strands[1]
    target_rstrand_models, rstrand_feature_stacks = strands[2]

    if output !== nothing
        if typeof(output) == String
            fname = output::String
            if isdir(fname)
                fname = joinpath(fname, "$(target_id).sff")
            end
        else
            fname = output # IOBuffer, IOStream
        end
    else
        fname = "$(target_id).sff"
    end

    # find inverted repeat if any
    f_aligned_blocks, r_aligned_blocks = alignLoops(targetloopf, target_saf, target_raf,
                                                    targetloopr, target_sar, target_rar)

    # sort blocks by length
    f_aligned_blocks = sort(f_aligned_blocks, by=last, rev=true)
    
    ir = if length(f_aligned_blocks) > 0 && f_aligned_blocks[1][3] >= 1000
        f_aligned_blocks[1]
    else
        nothing
    end
    
    writeSFF(fname, target_id, target_fstrand_models, target_rstrand_models,
        reference.gene_exons, fstrand_feature_stacks, rstrand_feature_stacks,
        targetloopf, targetloopr, ir)

    @info "[$target_id] Overall: $(ns(time_ns() - t1))"
    return fname, target_id

end

function annotate_one(reference::Reference, fasta::Union{String,IO})
    annotate_one(reference, fasta, IOBuffer())
end

function annotate(refsdir::String, templates::String, fa_files::Vector{String}, output::MayBeString)

    reference = readReferences(refsdir, templates)

    for infile in fa_files
        annotate_one(reference, infile, output)
    end
end

#### these are only used by chloe_distributed ####
function annotate_one_task(fasta::MayBeString, output::MayBeIO, task_id::MayBeString)
    annotation_local_storage(TASK_KEY, task_id)
    try
        # the global REFERENCE should have been
        # sent to the worker process by main process
        @debug "using $(Main.REFERENCE)"
        annotate_one(Main.REFERENCE::Reference, fasta, output)
    finally
        annotation_local_storage(TASK_KEY, nothing)
    end
end


function annotate_one_task(fasta::Union{String,IO}, task_id::MayBeString)
    annotation_local_storage(TASK_KEY, task_id)
    try
        # the global REFERENCE should have been
        # sent to the worker process by main process
        @debug "using $(Main.REFERENCE)"
        annotate_one(Main.REFERENCE::Reference, fasta, IOBuffer())
    finally
        annotation_local_storage(TASK_KEY, nothing)
    end
end
