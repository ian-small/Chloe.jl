include("../Utilities.jl")
include("../SuffixArrays.jl")
include("../Alignments3.jl")

target_id, target_seqf = readFasta(ARGS[1]) # "/Users/ian/github/pyrimid/AP000423.fa")
target_length = length(target_seqf)
targetloopf = target_seqf * target_seqf[1:end - 1]
target_saf = makeSuffixArray(targetloopf, true)
target_raf = makeSuffixArrayRanksArray(target_saf)

target_seqr = revComp(target_seqf)
targetloopr = target_seqr * target_seqr[1:end - 1]
target_sar = makeSuffixArray(targetloopr, true)
target_rar = makeSuffixArrayRanksArray(target_sar)

f_aligned_blocks, r_aligned_blocks = alignLoops(targetloopf, target_saf, target_raf, targetloopr, target_sar, target_rar)

# sort blocks by length
f_aligned_blocks = sort(f_aligned_blocks, by=last, rev=true)
ir = f_aligned_blocks[1]

if ir[3] >= 1000
    println("IR/1/repeat_region/1\t+\t" * string(ir[1]) * "\t" * string(ir[3]) * "\t0\t0\t0\t0\t")
    println("IR/2/repeat_region/1\t-\t" * string(ir[2]) * "\t" * string(ir[3]) * "\t0\t0\t0\t0\t")
end
