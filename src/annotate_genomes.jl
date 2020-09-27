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

datasize(m::MMappedString) = length(m.ptr)

struct SingleReference
    refsrc::String
    refloops::FwdRev{MappedPtrString}
    refSAs::FwdRev{SuffixArray}
    refRAs::FwdRev{SuffixArray}
    ref_features::Union{Nothing,FwdRev{FeatureArray}}
    target_length::Int32
    ismmapped::Bool
    forward_only::Bool
end
datasize(r::SingleReference) = begin
    (sizeof(Reference)
    + datasize(r.refsrc)
    + datasize(r.refloops)
    + datasize(r.refSAs)
    + datasize(r.refRAs)
    + datasize(r.ref_features)
    )
end
struct StrandSingleReference
    refsrc::String
    refloops::MappedPtrString
    refSAs::SuffixArray
    refRAs::SuffixArray
    ref_features::Union{Nothing,FeatureArray}
    target_length::Int32
    strand::Char
    ismmapped::Bool
end
datasize(r::StrandSingleReference) = begin
    (sizeof(StrandReference)
    + datasize(r.refsrc)
    + datasize(r.refloops)
    + datasize(r.refSAs)
    + datasize(r.refRAs)
    + datasize(r.ref_features)
    )
end
function stranded(strand::Char, r::SingleReference)::StrandSingleReference
    fa = r.ref_features
    if strand == '+'
        StrandSingleReference(r.refsrc, r.refloops.forward, r.refSAs.forward,
        r.refRAs.forward, fa === nothing ? nothing : fa.forward,
        r.target_length, strand, r.ismmapped)
    else
        StrandSingleReference(r.refsrc, r.refloops.reverse, r.refSAs.reverse,
        r.refRAs.reverse, fa === nothing ? nothing : fa.reverse, 
        r.target_length, strand, r.ismmapped)
    end
end

struct Reference
    references::Dict{String,SingleReference}
    # from .tsv file
    feature_templates::Dict{String,FeatureTemplate}
    gene_exons::Dict{String,Int32}

    # referenceOrganisms::Dict{String,String}
end
datasize(r::Reference) = begin
    (datasize(Reference) 
        + datasize(r.feature_templates) 
        + datasize(r.references) 
        + datasize(r.gene_exons)
    )
end

# stops the REPL printing the entire 7MB of sequences!
function Base.show(io::IO, reference::Reference)

    total = sum(datasize(s.second) for s in reference.references)
    bp = sum(datasize(s.second.refloops) for s in reference.references)
    # forward_only = sum(s.second.forward_only for s in reference.references)
    
    t1 = "#templates=$(reference.feature_templates |> length)"
    t2 = "#gene_exons=$(reference.gene_exons |> length)[$(reference.gene_exons |> values |> sum)]"
    t3 = "#refs=$(reference.references |> length)[$(human(bp)) bp]"
    # f =  "[forward#$(forward_only)]"
    print(io, "Reference: $t1, $t2, $t3, total=$(human(total))")

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
 
    if !isdir(refsdir)
        msg = "Reference directory $(refsdir) is not a directory!"
        @error msg
        throw(ArgumentError(msg))
    end

    for f in [template, joinpath(refsdir, "ReferenceOrganisms.json")]
        if !isfile(f)
            msg = "missing file: $f"
            @error msg
            throw(ArgumentError(msg))
        end
    end

    files = findall(x -> endswith(x, r"\.(gwsas|mmap)"), readdir(refsdir))
    if length(files) == 0
        msg = "please run `julia chloe.jl mmap $(refsdir)/*.{fa,fasta}`"
        @error msg
        throw(ArgumentError(msg))
    end
    
end

function readSingleReference(refsdir::String, ref::Pair{String,String}, 
    sff_files::AbstractArray{String};
    verbose::Bool=false,
    forward_only::Bool=false)
    gwsas = joinpath(refsdir, ref.first * ".gwsas")
    mmap = joinpath(refsdir, ref.first * ".mmap")
    ismmapped = false
    if isfile(mmap)
        refSAs, refRAs, refloops = read_mmap_suffix(mmap; forward_only=forward_only)
        ismmapped = true
        if verbose
            @info "found mmap file for: $(ref.first)"
        end
    elseif isfile(gwsas)
        refSAs, refRAs, refloops = read_gwsas(gwsas, ref.first; forward_only=forward_only)
        if verbose
            @info "found gwsas file for: $(ref.first)"
        end
    else
        error("no gwsas or mmap file found for \"$(ref.first)\". run `julia chloe.jl mmap --help`")
    end
    idx = findfirst(x -> startswith(x, ref.second), sff_files)
    if idx === nothing
        idx = findfirst(x -> startswith(x, ref.first), sff_files)
    end
    if idx === nothing
        # TODO check for fasta files
        error("no sff file for $(ref.first) -> $(ref.second)")
    end
    ref_features = readFeatures(joinpath(refsdir, sff_files[idx]))
    target_length = length(refSAs.forward)

    SingleReference(ref.first, refloops, refSAs, refRAs, ref_features, target_length, ismmapped, forward_only)
end

"""
    readReferences(reference_dir, template_file_tsv)

creates a Reference object to feed into `annotate_one`
"""
function readReferences(refsdir::String, templates::String;
    verbose::Bool=true, 
    forward_only::Bool=false)::Reference

    if !isdir(refsdir)
        error("$(refsdir) is not a directory")
    end
    path = joinpath(refsdir, "ReferenceOrganisms.json")
    if !isfile(path)
        error("readReferences: require \"$(path)\" JSON file")
    end
    ReferenceOrganisms = open(path) do f
        JSON.parse(f, dicttype=Dict{String,String})
    end

    num_refs = length(ReferenceOrganisms)
    refsdict = Dict{String,SingleReference}()

    files = readdir(refsdir)
    idx = findall(x -> endswith(x, ".sff"), files)
    if isempty(idx)
        error("No sff files found in $(refsdir)!")
    end
    
    sff_files = files[idx]
    mmaps = 0

    for ref in ReferenceOrganisms
        refsdict[ref.first] = r = readSingleReference(refsdir, ref, sff_files; verbose=verbose, forward_only=forward_only)
        if r.ismmapped 
            mmaps += 1
        end
    end
    if mmaps != length(refsdict)
        @info "better to use memory mapped files use: `julia chloe.jl mmap $(refsdir)/*.{fa,fasta}`"
    end
    feature_templates, gene_exons = readTemplates(templates)
    ret = Reference(refsdict, feature_templates, gene_exons)

    @info ret
    GC.gc() # cleanup
    return ret
end


MayBeString = Union{Nothing,String}
Strand = Tuple{AAFeature,DFeatureStack}
AAlignedBlocks = Vector{FwdRev{AlignedBlocks}}

function flatten(vanno::Vector{Vector{Annotation}})::Vector{Annotation}
    ret = Vector{Annotation}(undef, sum(length(v) for v in vanno))
    i = 1
    @inbounds for v in vanno
        for a in v
            ret[i] = a
            i += 1
        end
    end
    ret
end

function do_annotations(target_id::String, strand::Char, idxmap::Dict{Int,SingleReference}, blocks_aligned_to_target::AAlignedBlocks)
    # this takes about 4secs!
    function do_one(refsrc, ref_features, blocks)
        annotations = Vector{Annotation}()
        st = time_ns()
        push!(annotations, findOverlaps(ref_features.forward, blocks.forward)...)
        push!(annotations, findOverlaps(ref_features.reverse, blocks.reverse)...)
        @info "[$(target_id)]$(strand) $(refsrc)± overlaps $(length(annotations)): $(elapsed(st))"
        return annotations
    end
    tgt = Vector{Vector{Annotation}}(undef, length(idxmap))

    Threads.@threads for i in collect(idxmap)
        blocks = blocks_aligned_to_target[i.first]
        ref_features = i.second.ref_features
        tgt[i.first] = do_one(i.second.refsrc, ref_features, blocks)

    end
    # annotations = collect(Iterators.flatten(tgt))
    annotations = flatten(tgt)
    sort!(annotations, by=x -> x.path)
    annotations
end
function do_strand(target_id::String, start_ns::UInt64, target_length::Int32,
    idxmap::Dict{Int,SingleReference}, coverages::Dict{String,Float32},
    strand::Char, blocks_aligned_to_target::AAlignedBlocks,
    targetloop::DNAString, feature_templates::Dict{String,FeatureTemplate})::Strand

    annotations = do_annotations(target_id, strand, idxmap, blocks_aligned_to_target)

    t4 = time_ns()
    @info "[$target_id]$strand overlapping ref annotations ($(length(annotations))) $(human(datasize(annotations))): $(ns(t4 - start_ns))"

    # strand_feature_stacks is basically grouped by annotations.path
    strand_feature_stacks, shadow = fillFeatureStack(target_length, annotations, feature_templates)

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
            @debug "[$target_id]$strand below threshold: $(stack.path): depth=$(@sprintf "%.3f" depth) coverage=$( @sprintf "%.3f" coverage)"
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
    target_strand_models = refineGeneModels!(target_strand_models, target_length, 
                                            targetloop, annotations,
                                             path_to_stack)

    t8 = time_ns()
    @info "[$target_id]$strand refining gene models: $(ns(t8 - t7))"
    return target_strand_models, path_to_stack
end

# MayBeIO: write to file (String), IO buffer or create filename based on fasta filename
# Union{IO,String}: read fasta from IO buffer or a file (String)
MayBeIO = Union{String,IO,Nothing}

function createTargetReference(target_id::String, target_seqf::DNAString)::SingleReference
    target_seqr = revComp(target_seqf)
    target_length = Int32(length(target_seqf))
    
    targetloopf = target_seqf * target_seqf[1:end - 1]
    target_saf = makeSuffixArray(targetloopf, true)
    target_raf = makeSuffixArrayRanksArray(target_saf)

    targetloopr = target_seqr * target_seqr[1:end - 1]
    target_sar = makeSuffixArray(targetloopr, true)
    target_rar = makeSuffixArrayRanksArray(target_sar)
    
    targetloopf = MappedPtrString(targetloopf)
    targetloopr = MappedPtrString(targetloopr)
    SingleReference(target_id,
        FwdRev(targetloopf, targetloopr),
        FwdRev(target_saf, target_sar),
        FwdRev(target_raf, target_rar),
        nothing, # no .sff file
        target_length,
        false, false)
end

function alignit(tgt::SingleReference, ref::SingleReference)
    start = time_ns()
    src_id = ref.refsrc

    ff, fr = alignLoops(src_id, ref.refloops.forward, ref.refSAs.forward, ref.refRAs.forward, 
                                tgt.refloops.forward, tgt.refSAs.forward, tgt.refRAs.forward) 
    if !ref.forward_only
        rf, rr = alignLoops(src_id, ref.refloops.reverse, ref.refSAs.reverse, ref.refRAs.reverse, 
                                    tgt.refloops.forward, tgt.refSAs.forward, tgt.refRAs.forward)
    else
        rr, rf = alignLoops(src_id, ref.refloops.forward, ref.refSAs.forward, ref.refRAs.forward, 
                                    tgt.refloops.reverse, tgt.refSAs.reverse, tgt.refRAs.reverse)
    end

    # blocks_aligned_to_targetf[refcount] = FwdRev(ff, rf)
    # blocks_aligned_to_targetr[refcount] = FwdRev(rr, fr)

    @info "[$(tgt.refsrc)]± aligned $(src_id) ($(length(ff)),$(length(rf))) $(elapsed(start))"
    # note cross ...
    FwdRev(ff, rf), FwdRev(rr, fr)
end

function inverted_repeat(target::SingleReference)::AlignedBlock
    f_aligned_blocks, _ = alignLoops(target.refsrc,
                                                    target.refloops.forward, target.refSAs.forward, target.refRAs.forward,
                                                    target.refloops.reverse, target.refSAs.reverse, target.refRAs.reverse;
                                                    rev=false
                                        )

    # sort blocks by length
    f_aligned_blocks = sort!(f_aligned_blocks, by=last, rev=true)
    
    ir = if length(f_aligned_blocks) > 0
        f_aligned_blocks[1]
    else
        AlignedBlock(0, 0, 0)
    end
    ir
end

WatsonCrick = Tuple{Vector{Vector{SFF}},DFeatureStack}
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
function annotate_one(reference::Reference, target_id::String, target_seqf::String, output::MayBeIO)

    t1 = time_ns()

    target_length = Int32(length(target_seqf))
    n = count(r"[XN]", target_seqf)
    r = n / length(target_seqf)
    if r > .9
        error("sequence [$(target_id)] too vague: $(@sprintf "%.1f" r * 100)%  either X or N")
    end
    
    @info "[$target_id] seq length: $(target_length)bp"
    
    # target_seqr = revComp(target_seqf)
    
    # targetloopf = target_seqf * target_seqf[1:end - 1]
    # target_saf = makeSuffixArray(targetloopf, true)
    # target_raf = makeSuffixArrayRanksArray(target_saf)

    # targetloopr = target_seqr * target_seqr[1:end - 1]
    # target_sar = makeSuffixArray(targetloopr, true)
    # target_rar = makeSuffixArrayRanksArray(target_sar)
    
    # targetloopf = MappedPtrString(targetloopf)
    # targetloopr = MappedPtrString(targetloopr)

    target = createTargetReference(target_id, target_seqf)

    t2 = time_ns()

    @info "[$target_id] made suffix arrays $(human(datasize(target))): $(ns(t2 - t1))"

    num_refs = length(reference.references)

    blocks_aligned_to_targetf = AAlignedBlocks(undef, num_refs)
    blocks_aligned_to_targetr = AAlignedBlocks(undef, num_refs)
    idxmap = Dict{Int,SingleReference}()
    # map index into blocks_aligned_to_target{r,f} to SingleReference
    for (i, r) in enumerate(reference.references)
        idxmap[i] = r.second
    end


    Threads.@threads for i in collect(idxmap)
        fwd, rev = alignit(target, i.second)
        blocks_aligned_to_targetf[i.first] = fwd
        blocks_aligned_to_targetr[i.first] = rev
    end

    t3 = time_ns()
    
    @info "[$target_id] aligned: ($(num_refs)) $(human(datasize(blocks_aligned_to_targetf) + datasize(blocks_aligned_to_targetr))) $(ns(t3 - t2))" 

    coverages = Dict{String,Float32}()
    for key in idxmap
        a = blocks_aligned_to_targetf[key.first]
        coverage = blockCoverage(a.forward) + blockCoverage(a.reverse)
        coverages[key.second.refsrc] = coverage / (target.target_length * 2)
    end
    @debug "[$target_id] coverages:" coverages

    strands = Vector{WatsonCrick}(undef, 2)

    function watson()
        models, stacks = do_strand(target_id, t3, target.target_length, idxmap, coverages,
            '+', blocks_aligned_to_targetf, target.refloops.forward, reference.feature_templates)

        sffs_fwd = [toSFF(model, target.refloops.forward, stacks) 
            for model in filter(m -> !isempty(m), models)]
        strands[1] = (sffs_fwd, stacks)
        
    end
    function crick()
        models, stacks = do_strand(target_id, t3, target.target_length, idxmap, coverages,
            '-', blocks_aligned_to_targetr, target.refloops.reverse, reference.feature_templates)
        sffs_rev = [toSFF(model, target.refloops.reverse, stacks) 
            for model in filter(m -> !isempty(m), models)]
        strands[2] = (sffs_rev, stacks)
    end

    irb = Vector{AlignedBlock}(undef, 1)

    function bginverted_repeat()
        # find inverted repeat if any
        irb[1] = inverted_repeat(target)
    end

    Threads.@threads for worker in [watson, crick, bginverted_repeat]
        worker()
    end

    sffs_fwd, fstrand_feature_stacks = strands[1]
    sffs_rev, rstrand_feature_stacks = strands[2]
    ir = irb[1]
   
    
    if ir[3] >= 1000
        @info "[$target_id] inverted repeat $(ir[3])"
    else
        ir = nothing
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
    writeSFF(fname, target_id, target.target_length, reference.gene_exons,
             FwdRev(sffs_fwd, sffs_rev), ir)

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
    target_id, target_seqf = readFasta(fasta)
    annotate_one(reference, target_id, target_seqf, IOBuffer())
end

function annotate_all(reference::Reference, fasta::Union{String,IO})
    for (target_id, target_seqf) in iterFasta(fasta)
        annotate_one(reference, target_id, target_seqf, nothing)
    end
end

function annotate(refsdir::String, templates::String, fa_files::Vector{String}, output::MayBeString;
    verbose::Bool=true, forward_only::Bool=false)

    reference = readReferences(refsdir, templates; verbose=verbose, forward_only=forward_only)

    for infile in fa_files
        for (seq_id, seqf) in iterFasta(infile)
            annotate_one(reference, seq_id, seqf, output)
        end
    end
end

end # module
