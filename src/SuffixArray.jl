
module SuffixArray
export create_mmaps, writesuffixarray
import Mmap
import ..Annotator: revComp, readFasta, makeSuffixArray, makeSuffixArrayRanksArray, GenomeWithSAs, writeGenomeWithSAs

function writesuffixarray(;fasta_files=String[], directory=Union{String,Nothing} = nothing)
    for infile in fasta_files
        seqid, seqf = readFasta(infile)
        seqr = revComp(seqf)

        saf = makeSuffixArray(seqf * seqf[1:end - 1], true)
        sar = makeSuffixArray(seqr * seqr[1:end - 1], true)

        @assert length(saf) === length(sar)

        gwsas = GenomeWithSAs(seqid, seqf, saf, sar)
        if directory !== nothing
            filename = joinpath(directory, gwsas.id * ".gwsas")
        else
            d = splitpath(infile)
            filename = joinpath(d[1:end - 1]..., "$(gwsas.id).gwsas")

        end
        @info "writing suffix array: $filename"
        writeGenomeWithSAs(filename, gwsas)
    end
end

function one_mmap_from_gwsas(infile::String, forward_only::Bool=false)

    _, f = splitpath(infile)
    seqid = splitext(f)[1]
    gwsas = readGenomeWithSAs(infile, seqid)

    fwd = gwsas.sequence
    fwd = uppercase(fwd)
    fwds = fwd * fwd[1:end - 1]
    ufwd = Vector{UInt8}(fwds)
    
    function dorev()
        rev = revComp(fwd)
        revs = rev * rev[1:end - 1]
        Vector{UInt8}(revs)
    end
    
    n = length(gwsas.forwardSA)
    d = splitpath(infile)
    outfile = joinpath(d[1:end - 1]..., "$(seqid).mmap")
    f = forward_only ? "[forward-only]" : ""
    @info "writing mmap array$f from gwsas: $outfile"
    open(outfile, "w") do f
        write(f, gwsas.forwardSA)
        !forward_only && write(f, gwsas.reverseSA)
        write(f, makeSuffixArrayRanksArray(gwsas.forwardSA))
        !forward_only && write(f, makeSuffixArrayRanksArray(gwsas.reverseSA))
        write(f, ufwd)
        !forward_only && write(f, dorev())
    end
    @assert filesize(outfile) === (!forward_only ? 20 * n - 2 : 10 * n - 1)

end


function one_mmap_from_fasta(infile::String, forward_only::Bool=false)
    seqid, fwd = readFasta(infile)
    fwd = uppercase(fwd)

    fwds = fwd * fwd[1:end - 1]
    ufwd = Vector{UInt8}(fwds)

    saf = makeSuffixArray(fwds, true)

    if !forward_only
        rev = revComp(fwd)
        revs = rev * rev[1:end - 1]
    else
        rev = revs = ""
    end
    sar = makeSuffixArray(revs, true)
    urev = Vector{UInt8}(revs)

    
    n = length(saf)
    d = splitpath(infile)
    outfile = joinpath(d[1:end - 1]..., "$(seqid).mmap")
    f = forward_only ? "[forward-only]" : ""
    @info "writing mmap array$f from fasta to: $outfile"
    open(outfile, "w") do f
        write(f, saf)
        !forward_only && write(f, sar)
        write(f, makeSuffixArrayRanksArray(saf))
        !forward_only && write(f, makeSuffixArrayRanksArray(sar))
        write(f, ufwd)
        !forward_only && write(f, urev)
    end
    @assert filesize(outfile) === (!forward_only ? 20 * n - 2 : 10 * n - 1)

end

    function create_mmaps(;gwsas_fasta=String[], forward_only::Bool=false)
    for infile in gwsas_fasta
        if endswith(infile, ".gwsas")
            one_mmap_from_gwsas(infile, forward_only)
        elseif endswith(infile, r"\.(fasta|fa)")
            one_mmap_from_fasta(infile, forward_only)
        else
            @warn "unknown file type $infile"
        end
    end
end
end # module