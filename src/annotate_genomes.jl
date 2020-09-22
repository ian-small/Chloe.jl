module Annotator

export annotate, annotate_one, readReferences, Reference, MayBeIO, MayBeString
export MMappedString, ASCII

DNAString = AbstractString

import Base

include("UInt8Utf8.jl")
include("MMappedString.jl")
include("Utilities.jl")
include("SuffixArrays.jl")
include("Alignments3.jl")
include("Annotations.jl")

# import Dates: Time, Nanosecond

import Printf: @sprintf
import JSON
import Mmap
import Crayons: @crayon_str

import .MappedString: MMappedString, ASCII

const success = crayon"bold green"

# can be any of String, MMappedString, MappedPtrString{ASCII}
# If you have memory mapped files (see julia chloe.jl mmap *.fa)
# then use MMappedString{ASCII} since it will not read
# the data backing the memory mapping
MappedPtrString = MMappedString{ASCII}
# MappedPtrString = String

struct FwdRev{T}
    forward::T
    reverse::T
end

Base.:(==)(x::FwdRev{T}, y::FwdRev{T}) where T = x.forward == y.forward && x.reverse == y.reverse

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
    # feature_templates::Vector{FeatureTemplate}
    feature_templates::Dict{String,FeatureTemplate}
    gene_exons::Dict{String,Int32}
    forward_only::Bool
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
    if num == 0
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
    t3 = "#refs=$(reference.refloops |> length)[$(human(bp)) bp]"
    f = reference.forward_only ? "[forward]" : ""
    print(io, "Reference$(f): $t1, $t2, $t3, suffix=$(human(sabytes)), total=$(human(sabytes + bp))")

end

function read_mmap_suffix(filename::String; forward_only::Bool=false)
    size = filesize(filename)
    if isodd(size)
        if !forward_only
            error("memory mapped file \"$(filename)\" has been generated with forward sequences only (run with --forward-only)!")
        end

        return read_mmap_suffix_forward_only(filename)
    end
    # may as well mmap everything even thou we don't need it....

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
function read_mmap_suffix_forward_only(filename::String)
    size = filesize(filename)
    nelem = floor(Int, (size + 1) / 10)

    nbytes = 4 * nelem
    sbytes = 2 * nelem - 1
    zsa = SuffixArray()

    f = open(filename)
    # zero offset???
    saf = Mmap.mmap(f, SuffixArray, nelem, 0 + 0 * nbytes)
    # sar = Mmap.mmap(f, SuffixArray, nelem, 0 + 1 * nbytes)
    raf = Mmap.mmap(f, SuffixArray, nelem, 0 + 1 * nbytes)
    # rar = Mmap.mmap(f, SuffixArray, nelem, 0 + 3 * nbytes)
    sf =  Mmap.mmap(f, Vector{UInt8}, sbytes, 0 + 2 * nbytes)
    # sr  = Mmap.mmap(f, Vector{UInt8}, sbytes, 0 + 4 * nbytes + sbytes)
    msf, msr = MappedPtrString(sf), MappedPtrString(Vector{UInt8}())
    return FwdRev(saf, zsa), FwdRev(raf, zsa), FwdRev(msf, msr)
end

function read_gwsas(gwsas::String, id::String; forward_only::Bool=false)
    zsa = SuffixArray()   
    function FwdRevSA(gwsas::GenomeWithSAs)
        FwdRev(gwsas.forwardSA, forward_only ? zsa : gwsas.reverseSA)
    end
    function FwdRevRA(gwsas::GenomeWithSAs)
        FwdRev(makeSuffixArrayRanksArray(gwsas.forwardSA), forward_only ? zsa : makeSuffixArrayRanksArray(gwsas.reverseSA))
    end
    function FwdRevMap(fwd::String, rev::String)
        # makes a copy unfortuately...
        FwdRev(MappedPtrString(fwd), MappedPtrString(forward_only ? "" : rev))
    end
         
    refgwsas = readGenomeWithSAs(gwsas, id)
    fwd = refgwsas.sequence
    rev = revComp(refgwasas.sequence)
    
    fwd = fwd * fwd[1:end - 1]
    rev = rev * rev[1:end - 1]


    FwdRevSA(refgwsas), FwdRevRA(refgwsas), FwdRevMap(fwd, rev)
end

function verify_refs(refsdir, template)
    if ~isfile(template) || ~isdir(refsdir) || ~isfile(joinpath(refsdir, "ReferenceOrganisms.json"))
        msg = "template: $(template) or refsdir: $(refsdir) is incorrect!"
        @error msg
        throw(ArgumentError(msg))
    end
    files = findall(x -> endswith(x, r"\.(gwsas|mmap)"), readdir(refsdir))
    if length(files) == 0
        msg = "please run `julia chloe.jl mmap $(refsdir)/*.fa`"
        @error msg
        throw(ArgumentError(msg))
    end
    
end

"""
    readReferences(reference_dir, template_file_tsv)

creates a Reference object to feed into `annotate_one`
"""
function readReferences(refsdir::String, templates::String; verbose::Bool=true, 
    forward_only::Bool=false)::Reference

    if !isdir(refsdir)
        error("$(refsdir) is not a directory")
    end
    path = joinpath(refsdir, "ReferenceOrganisms.json")
    if ~isfile(path)
        error("readReferences: require \"$(path)\" JSON file")
    end
    ReferenceOrganisms = open(path) do f
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
    mmaps = 0

    for (i, ref) in enumerate(ReferenceOrganisms)
        gwsas = joinpath(refsdir, ref.first * ".gwsas")
        mmap = joinpath(refsdir, ref.first * ".mmap")
        
        if isfile(mmap)
            refSAs[i], refRAs[i], refloops[i] = read_mmap_suffix(mmap; forward_only=forward_only)
            if verbose
                @info "found mmap file for: $(ref.first)"
            end
            mmaps += 1
        elseif isfile(gwsas)
            refSAs[i], refRAs[i], refloops[i] = read_gwsas(gwsas, ref.first; forward_only=forward_only)
            if verbose
                @info "found gwsas file for: $(ref.first)"
            end
        else
            error("no data file for $(ref.first)")
        end
        idx = findfirst(x -> startswith(x, ref.second), sff_files)
        if idx === nothing
            idx = findfirst(x -> startswith(x, ref.first))
        end
        if idx === nothing
            # TODO check for fasta files
            error("no sff file for $(ref.first) -> $(ref.second)")
        end
        feature_file = sff_files[idx]
        f_strand_features, r_strand_features = readFeatures(joinpath(refsdir, feature_file))

        ref_features[i] = FwdRev(f_strand_features, r_strand_features)
        refsrc[i] = ref.first

    end
    if mmaps != length(refsrc)
        @info "better to use memory mapped files use: `julia chloe.jl mmap $(refsdir)/*.fa`"
    end
    feature_templates, gene_exons = readTemplates(templates)
    ret = Reference(refsrc, refloops, refSAs, refRAs, ref_features, feature_templates, gene_exons, forward_only)

    @info ret
    GC.gc() # cleanup
    return ret
end


# const ns(td) = Time(Nanosecond(td))
const ns(td) = @sprintf("%.3fs", td / 1e9)

MayBeString = Union{Nothing,String}
Strand = Tuple{AAFeature,DFeatureStack}
AAlignedBlocks = Vector{FwdRev{AlignedBlocks}}

function do_strand(target_id::String, start_ns::UInt64, target_length::Int32,
    reference::Reference, coverages::Dict{String,Float32},
    strand::Char, blocks_aligned_to_target::AAlignedBlocks,
    targetloop::DNAString)::Strand

    annotations = Vector{Annotation}()
    for (ref_feature_array, blocks) in zip(reference.ref_features, blocks_aligned_to_target)
        push!(annotations, findOverlaps(ref_feature_array.forward, blocks.forward)...)
        push!(annotations, findOverlaps(ref_feature_array.reverse, blocks.reverse)...)
    end
    sort!(annotations, by=x -> x.path)

    t4 = time_ns()
    @info "[$target_id]$strand overlapping ref annotations ($(length(annotations))) $(human(datasize(annotations))): $(ns(t4 - start_ns))"

    # strand_feature_stacks is basically grouped by annotations.path
    strand_feature_stacks, shadow = fillFeatureStack(target_length, annotations, reference.feature_templates)

    t5 = time_ns()
    @info "[$target_id]$strand ref features stacks ($(length(strand_feature_stacks))) $(human(datasize(strand_feature_stacks))): $(ns(t5 - t4))"

    features = Vector{Feature}()
    # sizehint!(fa, length(strand_feature_stacks))
    path_to_stack = Dict{String,FeatureStack}()
    # so... features will be ordered by annotations.path
    for stack in strand_feature_stacks
        left_border, length = alignTemplateToStack(stack, shadow)
        left_border == 0 && continue
        depth, coverage = getDepthAndCoverage(stack, left_border, length)
        if ((depth >= stack.template.threshold_counts) && (coverage >= stack.template.threshold_coverage))
            push!(features, Feature(stack.path, left_border, length, 0))
            if haskey(path_to_stack, stack.path)
                @error "duplicate feature path: $(stack.path)"
            end
            path_to_stack[stack.path] = stack
        else
            @debug "[$target_id]$strand Below threshold: $(stack.path)"
        end
    end

    t6 = time_ns()
    @info "[$target_id]$strand aligning templates ($(length(features))): $(ns(t6 - t5))"

    for feature in features
        refineMatchBoundariesByOffsets!(feature, annotations, target_length, coverages)
    end

    t7 = time_ns()
    @info "[$target_id]$strand refining match boundaries: $(ns(t7 - t6))"

    # group by feature name on **ordered** features getFeatureName()
    target_strand_models = groupFeaturesIntoGeneModels(features)
    # this toys with the feature start, phase etc....
    target_strand_models = refineGeneModels!(target_strand_models, target_length, targetloop, annotations,
                                             path_to_stack)

    t8 = time_ns()
    @info "[$target_id]$strand refining gene models: $(ns(t8 - t7))"
    return target_strand_models, path_to_stack
end

# MayBeIO: write to file (String), IO buffer or create filename based on fasta filename
# Union{IO,String}: read fasta from IO buffer or a file (String)
MayBeIO = Union{String,IO,Nothing}

"""
    annotate_one(references::Reference, fasta_file::Union{String,IO} [,output_sff_file])

Annotate a single fasta file containting a *single* circular
DNA entry

writes an .sff file to `output_sff_file` or uses the sequence id in the
fasta file to write `{target_id}.sff` in the current directory.

If `output_sff_file` is a *directory* write `{target_id}.sff` into that
directory.

returns a 2-tuple: (ultimate sff output filename, sequence id)

If `output_sff_file` is an IOBuffer then that buffer will be returned
with the annotation within it

`reference` are the reference annotations (see `readReferences`)
"""
function annotate_one(reference::Reference, fasta::Union{String,IO}, output::MayBeIO)

    num_refs = length(reference.refsrc)
    t1 = time_ns()

    target_id, target_seqf = readFasta(fasta)
    target_length = Int32(length(target_seqf))
    
    @info "[$target_id] seq length: $(target_length)bp"
    
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

    @info "[$target_id] made suffix arrays $(human(datasize(target_saf) + datasize(target_raf) + datasize(target_sar) + datasize(target_rar))): $(ns(t2 - t1))"

    blocks_aligned_to_targetf = AAlignedBlocks(undef, num_refs)
    blocks_aligned_to_targetr = AAlignedBlocks(undef, num_refs)

    function alignit(refcount::Int)
        start = time_ns()
        refloop, refSA, refRA = reference.refloops[refcount], reference.refSAs[refcount], reference.refRAs[refcount]

        ff, fr = alignLoops(refloop.forward, refSA.forward, refRA.forward, targetloopf, target_saf, target_raf) 
        if !reference.forward_only
            rf, rr = alignLoops(refloop.reverse, refSA.reverse, refRA.reverse, targetloopf, target_saf, target_raf)
        else
            rr, rf = alignLoops(refloop.forward, refSA.forward, refRA.forward, targetloopr, target_sar, target_rar)
        end
        # note cross ...
        blocks_aligned_to_targetf[refcount] = FwdRev(ff, rf)
        blocks_aligned_to_targetr[refcount] = FwdRev(rr, fr)
    
        @info "[$target_id]± aligned $(reference.refsrc[refcount]) ($(length(ff)),$(length(rf))) $(ns(time_ns() - start))"

    end

    Threads.@threads for refno in 1:length(reference.refloops)
        alignit(refno)
    end

    t3 = time_ns()
    
    @info "[$target_id] aligned: ($(length(reference.refloops))) $(human(datasize(blocks_aligned_to_targetf) + datasize(blocks_aligned_to_targetr))) $(ns(t3 - t2))" 

    coverages = Dict{String,Float32}()
    for (ref, a) in zip(reference.refsrc, blocks_aligned_to_targetf)
        coverage = blockCoverage(a.forward) + blockCoverage(a.reverse)
        coverages[ref] = coverage / target_length * 2
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
    
    if output !== nothing
        if typeof(output) == String
            fname = output::String
            if isdir(fname)
                fname = joinpath(fname, "$(target_id).sff")
            end
        else
            fname = output # IOBuffer
        end
    else
        fname = "$(target_id).sff"
    end
    writeSFF(fname, target_id, target_length, reference.gene_exons,
        target_fstrand_models, target_rstrand_models,
        fstrand_feature_stacks, rstrand_feature_stacks,
        targetloopf, targetloopr,
        ir)

    @info success("[$target_id] Overall: $(ns(time_ns() - t1))")
    return fname, target_id

end
"""
    annotate_one(reference::Reference, fasta::Union{String,IO})

Annotate a fasta file. Maybe a file name or an IOBuffer
returns a 2-tuple (sff annotation as an IOBuffer, sequence id)

`reference` are the reference annotations (see `readReferences`)
"""
function annotate_one(reference::Reference, fasta::Union{String,IO})
    annotate_one(reference, fasta, IOBuffer())
end

function annotate(refsdir::String, templates::String, fa_files::Vector{String}, output::MayBeString;
    verbose::Bool=true, forward_only::Bool=false)

    reference = readReferences(refsdir, templates; verbose=verbose, forward_only=forward_only)

    for infile in fa_files
        annotate_one(reference, infile, output)
    end
end

end # module
