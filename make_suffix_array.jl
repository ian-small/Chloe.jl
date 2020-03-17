include("Utilities.jl")
include("SuffixArrays.jl")

for infile in ARGS
    seqid, seqf = readFasta(infile)
    seqr = revComp(seqf)

    saf = makeSuffixArray(seqf*seqf[1:end-1],true)
    sar = makeSuffixArray(seqr*seqr[1:end-1],true)

    gwsas = GenomeWithSAs(seqid,seqf,saf,sar)
    filename = gwsas.id * ".gwsas"
    writeGenomeWithSAs(filename,gwsas)
end
