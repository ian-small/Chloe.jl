module Annotator

export annotate, annotate_one, readReferences, Reference, MayBeIO, MayBeString
export readSingleReference, inverted_repeat

import Base

include("Utilities.jl")
include("hash.jl")
include("AlignGenomes.jl")
include("orfs.jl")
include("Annotations.jl")

# import Dates: Time, Nanosecond

import Printf: @sprintf
import JSON
import Crayons: @crayon_str
using StatsBase
using Flux

const success = crayon"bold green"
const MINIMUM_UNASSIGNED_ORF_LENGTH = Int32(270)

struct SingleReference
    ref_id::String
    ref_seq::CircularSequence
    ref_features::Union{Nothing,FwdRev{FeatureArray}}
end

datasize(r::SingleReference) = begin
    (sizeof(Reference)
    + datasize(r.refsrc)
    + datasize(r.ref_features)
    )
end

function Base.show(io::IO, r::SingleReference)
    total = datasize(r)
    bp = datasize(length(r.refsrc))
    print(io, "SingleReference[$(r.refsrc)]: $(human(bp))bp, total=$(human(total))")
end

function Base.show(io::IO, r::FwdRev{FwdRev{AlignedBlocks}})
    total = datasize(r)
    ff = length(r.forward.forward)
    fr = length(r.forward.reverse)
    rf = length(r.reverse.forward)
    rr = length(r.reverse.reverse)
    print(io, "Chloë Alignment: [($ff,$fr),($rf,$rr)] total=$(human(total))")
end

function verify_refs(refsdir, template)
    # used by master process to check reference directory
    # *before* starting worker processes...

    # TODO: read json file and really check...
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
    files = readdir(refsdir)    
    sff = findall(x -> endswith(x, ".sff"), files)
    if length(sff) == 0
        msg = "no reference .sff files in $(refsdir)!"
        @error msg
        throw(ArgumentError(msg))
    end
end

function readSingleReference(refdir::String, refID::AbstractString)::SingleReference
    ref = FASTA.Record()
    if !isdir(refdir); refdir = dirname(refdir); end
    path = findfastafile(refdir, refID)
    if isnothing(path)
        msg = "unable to find $(refID) fasta file in $(refdir)!"
        @error msg
        throw(ArgumentError(msg))
    end
    reader = open(FASTA.Reader, path)
    read!(reader, ref)
    close(reader)
    ref_features = readFeatures(normpath(joinpath(refdir, refID * ".sff")))
    SingleReference(refID, CircularSequence(FASTA.sequence(ref)), ref_features)
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

function do_annotations(numrefs:: Int, target_id::String, strand::Char, refs::Vector{SingleReference}, blocks_aligned_to_target::AAlignedBlocks)

    function do_one(refsrc, ref_features, blocks)
        st = time_ns()
        annotations = findOverlaps(ref_features.forward, blocks.forward)
        annotations = vcat(annotations, findOverlaps(ref_features.reverse, blocks.reverse))
        @debug "[$(target_id)]$(strand) $(refsrc)± overlaps $(length(annotations)): $(elapsed(st))"
        return annotations
    end
    tgt = Vector{Vector{Annotation}}(undef, numrefs)

    Threads.@threads for i in 1:numrefs
        blocks = blocks_aligned_to_target[i]
        tgt[i] = do_one(refs[i].ref_seq, refs[i].ref_features, blocks)
    end

    # annotations = collect(Iterators.flatten(tgt))
    annotations = flatten(tgt)
    sort!(annotations, by = x -> x.path)
    annotations
end

function do_strand(numrefs::Int, target_id::String, target_seq::CircularSequence, refs::Vector{SingleReference}, coverages::Dict{String,Float32},
    strand::Char, blocks_aligned_to_target::AAlignedBlocks, feature_templates::Dict{String,FeatureTemplate})::Strand

    t4 = time_ns()
    annotations = do_annotations(numrefs, target_id, strand, refs, blocks_aligned_to_target)

    # strand_feature_stacks is basically grouped by annotations.path
    target_length = Int32(length(target_seq))
    strand_feature_stacks, shadow = fillFeatureStack(target_length, annotations, feature_templates)

    t5 = time_ns()
    @info "[$target_id]$strand built feature stacks ($(length(strand_feature_stacks))) $(human(datasize(strand_feature_stacks))): $(ns(t5 - t4))"

    #orfmap = threeframes(target_seq)
    features = Vector{Feature}()
    path_to_stack = Dict{String,FeatureStack}()
    # so... features will be ordered by annotations.path

    for stack in strand_feature_stacks
        left_border, length = alignTemplateToStack(stack, shadow)
        left_border == 0 && continue
        push!(features, Feature(stack.path, left_border, length, 0))
        if haskey(path_to_stack, stack.path)
            @error "duplicate feature path: $(stack.path)"
        end
        path_to_stack[stack.path] = stack
    end

    # t6 = time_ns()
    # @info "[$target_id]$strand transferring annotations ($(length(features))): $(ns(t6 - start_ns))"

    for feature in features
        refineMatchBoundariesByOffsets!(feature, annotations, target_length, coverages)
    end

    # t7 = time_ns()
    # @info "[$target_id]$strand refining match boundaries: $(ns(t7 - t6))"

    # group by feature name on **ordered** features getFeatureName()
    target_strand_models = groupFeaturesIntoGeneModels(features)
    # this toys with the feature start, phase etc....
    target_strand_models = refineGeneModels!(target_strand_models, target_seq, annotations, path_to_stack)

    #add any unassigned orfs
    uorfs = getallorfs(target_seq, strand, MINIMUM_UNASSIGNED_ORF_LENGTH)
    for uorf in uorfs
        #println(uorf)
        unassigned = true
        orf_frame = mod1(uorf.start, 3)
        for f in features
            #println(f)
            f_frame = mod1(f.start + f.phase, 3)
            if orf_frame == f_frame && rangesOverlap(uorf.start, uorf.length, f.start, f.length)
                unassigned = false
            end
        end
        if unassigned
            #count codons
            codonfrequencies = countcodons(uorf,target_seq)
            #predict with FluxCodingClassifier: 1 = coding, 2 = non-coding
            prediction = Flux.onecold(FluxCodingClassifier(codonfrequencies))
            if prediction == 1
                #predict with FluxGeneClassifier
                probs = FluxGeneClassifier(codonfrequencies)
                if maximum(probs) > 0.999 #sanity check to avoid daft predictions
                    prediction = Flux.onecold(probs)
                    #println(exon_indices[prediction],'\t', maximum(probs))
                    newtags = split(exon_indices[prediction], "/")
                    uorf._path_components[1] = "predicted_"*newtags[1]
                    uorf._path_components[4] = newtags[2]
                    uorf.path = annotationPath(uorf)
                end
                push!(target_strand_models, [uorf])
            end
        end
    end

    @info "[$target_id]$strand built gene models: $(elapsed(t5))"
    return target_strand_models, path_to_stack
end

# MayBeIO: write to file (String), IO buffer or create filename based on fasta filename
# Union{IO,String}: read fasta from IO buffer or a file (String)
MayBeIO = Union{String,IO,Nothing}

function align(ref::SingleReference, tgt_id::String, tgt_seq::CircularSequence, rev_tgt_seq::CircularSequence)::Tuple{FwdRev{FwdRev{AlignedBlocks}},Float32}
    start = time_ns()
    #align2seqs returns linked list
    #but subsequent functions expecting Vector
    ff::AlignedBlocks = ll2vector(align2seqs(ref.ref_seq.sequence, tgt_seq.sequence))
    rr = revCompBlocks(ff,length(ref.ref_seq.sequence), length(tgt_seq.sequence))
    
    fr::AlignedBlocks = ll2vector(align2seqs(ref.ref_seq.sequence,rev_tgt_seq.sequence))
    rf = revCompBlocks(fr,length(ref.ref_seq.sequence), length(tgt_seq.sequence))

    coverage = Float32(100 * target_coverage(ff, rf, length(tgt_seq)))
    #println(tgt_id,"\t",ref.ref_id,'\t', coverage)

    @info "[$(tgt_id)]± aligned $(ref.ref_id) coverage: $(@sprintf("%.2f", coverage))% $(elapsed(start))"
    # note cross ...
    (FwdRev(FwdRev(ff, rf), FwdRev(fr, rr)), coverage)
end

function inverted_repeat(target::CircularSequence, revtarget::CircularSequence)::AlignedBlock
    fr::AlignedBlocks = ll2vector(align2seqs(target.sequence,revtarget.sequence))
    # sort blocks by length
    fr = sort!(fr, by= b -> b.blocklength, rev=true)
    ir = length(fr) > 0 ? fr[1] : AlignedBlock(0, 0, 0)
    return ir
end

function avg_coverage(target::SingleReference, a::FwdRev{AlignedBlocks})
    coverage = blockCoverage(a.forward) + blockCoverage(a.reverse)
    coverage /= (target.target_length * 2)
    coverage
end

function avg_coverage(target::SingleReference, a::FwdRev{FwdRev{AlignedBlocks}})
    avg_coverage(target, a.forward)
end

Models = Vector{Vector{SFF}}

"""
    annotate_one(references::Reference, seq_id::String, seq::String, [,output_sff_file])

Annotate a single sequence containting a *single* circular
DNA entry

writes an .sff file to `output_sff_file` or uses the sequence id in the
fasta file to write `{seq_id}.sff` in the current directory.

If `output_sff_file` is a *directory* write `{seq_id}.sff` into that
directory.

returns a 2-tuple: (ultimate sff output filename, sequence id)

If `output_sff_file` is an IOBuffer then that buffer will be returned
with the annotation within it

`reference` are the reference annotations (see `readReferences`)
"""
function annotate_one(refsdir::String, numrefs::Int, refhashes::Union{Nothing,Dict{String,Vector{Int64}}}, templates_file::String, target_id::String,
     target_forward_strand::CircularSequence, target_reverse_strand::CircularSequence, output::MayBeIO)::Tuple{Union{String,IO},String}

    t1 = time_ns()

    fname = if output !== nothing
        if output isa String
            if isdir(output)
                joinpath(output, "$(target_id).sff")
            else
                output # filename
            end
        else
            output # IOBuffer
        end
    else
        "$(target_id).sff"
    end

    # sanity checks
    target_length = length(target_forward_strand)
    n = count(isambiguous, target_forward_strand.sequence)
    r = n / target_length
    if r > .01
        @warn("sequence [$(target_id)] contains too many ambiguous nucleotides: $(@sprintf "%.1f" r * 100)% ")
        return fname, target_id
    end
    n = count(isgap, target_forward_strand.sequence)
    if n > 0
        @warn("sequence [$(target_id)] contains gaps which will be removed before analysis")
        target_forward_strand = ungap!(target_forward_strand)
        target_reverse_strand = ungap!(target_reverse_strand)
        target_length = Int32(length(target_forward_strand))
    end
    
    @info "[$target_id] seq length: $(target_length)bp"

    # find best references
    if !isnothing(refhashes)
        hash = minhash(target_forward_strand.sequence, KMERSIZE, SKETCHSIZE)
        numrefs = min(numrefs, length(refhashes))
        refpicks = searchhashes(hash, refhashes)[1:numrefs]
    else
        numrefs = 1
        refpicks = [(split(basename(refsdir),'.')[1],0)]
    end

    t2 = time_ns()

    @info "[$target_id] picked $(numrefs) reference(s): $(ns(t2 - t1))"

    blocks_aligned_to_targetf = AAlignedBlocks(undef, numrefs)
    blocks_aligned_to_targetr = AAlignedBlocks(undef, numrefs)
    coverages = Dict{String,Float32}()
    refs = Vector{SingleReference}(undef, numrefs)

    Threads.@threads for i in 1:numrefs
        ref = readSingleReference(refsdir, refpicks[i][1])
        a = align(ref, target_id, target_forward_strand, target_reverse_strand)
        blocks_aligned_to_targetf[i] = a[1].forward
        blocks_aligned_to_targetr[i] = a[1].reverse
        lock(REENTRANT_LOCK)
        coverages[refpicks[i][1]] = a[2]
        unlock(REENTRANT_LOCK)
        refs[i] = ref
    end

    t3 = time_ns()
    @info "[$target_id] aligned: ($(numrefs)) $(human(datasize(blocks_aligned_to_targetf) + datasize(blocks_aligned_to_targetr))) mean coverage: $(geomean(values(coverages))) $(ns(t3 - t2))" 

    feature_templates = readTemplates(templates_file)

    function watson()
        models, stacks = do_strand(numrefs, target_id, target_forward_strand, refs, coverages,
            '+', blocks_aligned_to_targetf, feature_templates)

        [toSFF(model, target_forward_strand, stacks, numrefs) 
            for model in filter(m -> !isempty(m), models)]
    end

    function crick()
        models, stacks = do_strand(numrefs, target_id, target_reverse_strand, refs, coverages,
        '-', blocks_aligned_to_targetr, feature_templates)

        [toSFF(model, target_reverse_strand, stacks, numrefs) 
            for model in filter(m -> !isempty(m), models)]
    end

    # from https://discourse.julialang.org/t/threads-threads-to-return-results/47382

    sffs_fwd, sffs_rev, ir = fetch.((Threads.@spawn w()) for w in [watson, crick, () -> inverted_repeat(target_forward_strand, target_reverse_strand)]) 
    
    # sffs_fwd = watson()
    # sffs_rev = crick()
    # ir = inverted_repeat(target_seq, rev_target_seq)
   
    if ir.blocklength >= 1000
        @info "[$target_id] inverted repeat $(ir.blocklength)"
    else
        ir = nothing
    end
    
    writeSFF(fname, target_id, target_length, FwdRev(sffs_fwd, sffs_rev), ir)

    @info success("[$target_id] Overall: $(elapsed(t1))")
    return fname, target_id

end
"""
    annotate_one(reference::Reference, fasta::Union{String,IO})

Annotate a fasta file. Maybe a file name or an IOBuffer
returns a 2-tuple (sff annotation as an IOBuffer, sequence id)

`reference` are the reference annotations (see `readReferences`)
"""
# function annotate_one(reference::Reference, fasta::Union{String,IO})::Tuple{Union{String,IO},String}
#     target_id, target_seqf = readFasta(fasta)
#     annotate_one(reference, target_id, target_seqf, IOBuffer())
# end
# function annotate_one(reference::Reference, fasta::Union{String,IO}, output::MayBeIO)::Tuple{Union{String,IO},String}
#     target_id, target_seqf = readFasta(fasta)
#     annotate_one(reference, target_id, target_seqf, output)
# end

# function annotate_all(reference::Reference, fasta::Union{String,IO})
#     for (target_id, target_seqf) in iterFasta(fasta)
#         annotate_one(reference, target_id, target_seqf, nothing)
#     end
# end

function annotate(refsdir::String, numrefs::Int, hashfile::Union{String,Nothing}, templates::String, fa_files::Vector{String}, output::MayBeString)

    numrefs == 1 ? refhashes = nothing : refhashes = readminhashes(hashfile)

    for infile in fa_files
        reader = open(FASTA.Reader, infile)
        records = Vector{FASTA.Record}()
        for record in reader
            push!(records, record)
        end
        if isempty(records)
            @error "unable to read sequence from $infile"
        elseif length(records) > 2
            @error "$infile contains multiple sequences; Chloë expects a single sequence per file"
        end
        fseq = CircularSequence(FASTA.sequence(records[1]))
        if length(records) == 2
            rseq = CircularSequence(FASTA.sequence(records[2]))
        else
            rseq = reverse_complement(fseq)
        end
        if length(fseq) != length(rseq)
            @error "unequal lengths for forward and reverse strands"
        end
        annotate_one(refsdir, numrefs, refhashes, templates, FASTA.identifier(records[1]), fseq, rseq, output)
        close(reader)
    end
end

end # module

function allonall()
	Logging.disable_logging(Logging.Warn)
	files = filter(x->endswith(x, ".fa"), readdir("/Users/ian/github/chloe_references/", join = true))
	for file1 in files
		id1 = split(basename(file1),"/")[1]
		for file2 in files
			id2 = split(basename(file2),"/")[1]
			o = join(["/Users/ian/Data/cp_genomes/refseq/Spermatophyta/oneonone/",id1,"_",id2,".sff"])
			Annotator.annotate(file1, 1, nothing, "/Users/ian/github/chloe_references/empty_templates.tsv", [file2], o)
		end
	end
end