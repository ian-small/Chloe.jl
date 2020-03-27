# include("Utilities.jl")
# include("SuffixArrays.jl")

function writesuffixarray(;fasta_files = String[])
    for infile in fasta_files
        seqid, seqf = readFasta(infile)
        seqr = revComp(seqf)

        saf = makeSuffixArray(seqf * seqf[1:end - 1], true)
        sar = makeSuffixArray(seqr * seqr[1:end - 1], true)

        gwsas = GenomeWithSAs(seqid, seqf, saf, sar)
        filename = gwsas.id * ".gwsas"
        @info "writing suffix array: $(filename)"
        writeGenomeWithSAs(filename, gwsas)
    end
end
