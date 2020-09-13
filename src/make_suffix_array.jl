
import Mmap

function writesuffixarray(;fasta_files=String[], directory=Union{String,Nothing})
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
            filename = gwsas.id * ".gwsas"
        end
        @info "writing suffix array: $filename"
        writeGenomeWithSAs(filename, gwsas)
    end
end

function one_mmap_from_gwsas(infile::String)

    _, f = splitpath(infile)
    id = splitext(f)[1]
    gwsas = readGenomeWithSAs(infile, id)

    fwd = gwsas.sequence
    fwd = uppercase(fwd)
    rev = revComp(fwd)
    fwds = fwd * fwd[1:end - 1]
    revs = rev * rev[1:end - 1]

    ufwd = Vector{UInt8}(fwds)
    urev = Vector{UInt8}(revs)
    n = length(gwsas.forwardSA)
    d = splitpath(infile)
    outfile = joinpath(d[1:end - 1]..., "$(id).mmap")
    @info "writing mmap array from gwsas: $outfile"
    total = 4 * 4 * n  + 2 * ( 2 * n - 1) 
    open(outfile, "w") do f
        write(f, gwsas.forwardSA)
        write(f, gwsas.reverseSA)
        write(f, makeSuffixArrayRanksArray(gwsas.forwardSA))
        write(f, makeSuffixArrayRanksArray(gwsas.reverseSA))
        write(f, ufwd)
        write(f, urev)
    end
    @assert filesize(outfile) === 20 * n - 2

end


function one_mmap_from_fasta(infile::String)
    seqid, fwd = readFasta(infile)
    fwd = uppercase(fwd)
    rev = revComp(fwd)
    fwds = fwd * fwd[1:end - 1]
    revs = rev * rev[1:end - 1]
    saf = makeSuffixArray(fwds, true)
    sar = makeSuffixArray(revs, true)

    @assert length(saf) === length(sar)

    ufwd = Vector{UInt8}(fwds)
    urev = Vector{UInt8}(revs)
    n = length(saf)
    d = splitpath(infile)
    outfile = joinpath(d[1:end - 1]..., "$(seqid).mmap")
    @info "writing mmap array from fasta to: $outfile"
    total = 4 * 4 * n  + 2 * ( 2 * n - 1) 
    open(outfile, "w") do f
        write(f, saf)
        write(f, sar)
        write(f, makeSuffixArrayRanksArray(saf))
        write(f, makeSuffixArrayRanksArray(sar))
        write(f, ufwd)
        write(f, urev)
    end
    @assert filesize(outfile) === 20 * n - 2

end

function create_mmaps(;gwsas_fasta=String[])
    for infile in gwsas_fasta
        if endswith(infile, ".gwsas")
            one_mmap_from_gwsas(infile)
        elseif endswith(infile, ".fa")
            one_mmap_from_fasta(infile)
        else
            @warn "unknown file type $infile"
        end
    end
end