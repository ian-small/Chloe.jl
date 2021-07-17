using FASTX
include("sais.jl")
include("ChainedBlocks.jl")

global const MINIMUMFILLABLEGAP = 20
global const MAXIMUMMERGEABLEGAP = 400

function align(;query, target, output=Base.stdout)
    reader = open(FASTA.Reader, query)
    record = FASTA.Record()
    read!(reader, record)
    seq1 = FASTA.sequence(record)
    close(reader)
    for target_file in target
        reader = open(FASTA.Reader, target_file)
        record = FASTA.Record()
        read!(reader, record)
        seq2 = FASTA.sequence(record)
        alignment = Alignment(seq1, seq2)
        cseq1 = CircularSequence(seq1)
        emask = entropy_mask(cseq1, Int32(KMERSIZE))
        chain = align2seqs(cseq1, CircularSequence(seq2), emask)
        for link in chain
            block = link.data
            query_segment = LongDNASeq(collect(reinterpret(DNA, alignment[i]) for i in block.src_index:block.src_index + block.blocklength - 1))
            target_segment = LongDNASeq(collect(reinterpret(DNA, alignment[i]) for i in fulcrum(alignment) + block.tgt_index:fulcrum(alignment) + block.tgt_index + block.blocklength - 1))
            println("$(block.src_index)\t$query_segment\t$(block.src_index + block.blocklength - 1)")
            println("$(block.tgt_index)\t$target_segment\t$(block.tgt_index + block.blocklength - 1)\n")
        end
        close(reader)
    end
end

@inline function is_source(index::Int32, fulcrum::Int32)
    return index < fulcrum ? true : false
end

@inline function is_target(index::Int32, fulcrum::Int32)
    return index > fulcrum ? true : false
end

function align2seqs(seq1::CircularSequence, seq2::CircularSequence, mask=nothing)::BlockChain{AlignedBlock}
    # extend and join strings; currently assumes both sequence are circular
    # TO DO: generalise for either or both of seq1, seq2 as linear sequences
    alignment = Alignment(seq1.sequence, seq2.sequence)
    sa = sais(alignment, zeros(Int32, length(alignment)), 0, length(alignment), 16) .+= 1 # shift to 1-based julia indexes
    ra = invperm(sa)
    lcps = lcparray(alignment, sa, ra)
    lenseq1 = Int32(length(seq1))
    lenseq2 = Int32(length(seq2))

    src2tgt = blockchain(alignment, sa, ra, lcps, one(Int32):lenseq1, one(Int32):lenseq2, matchLengthThreshold(lenseq1, lenseq2), mask)
    src2tgt = circularise(src2tgt, lenseq1)

    #= println("links before gapfilling: $(src2tgt.links)")
    for link in src2tgt
        println(link)
    end =#

    Threads.@threads for gap in gaps(src2tgt)
        (srcgap, tgtgap) = contiguousblockgaps(gap[1], gap[2], lenseq1, lenseq2)
        #println(srcgap, " ", length(srcgap), "\t", tgtgap, " ", length(tgtgap))
        lengthsrcgap::Int32 = srcgap.stop - srcgap.start + 1
        lengthtgtgap::Int32 = tgtgap.stop - tgtgap.start + 1
        (lengthsrcgap < MINIMUMFILLABLEGAP || lengthtgtgap < MINIMUMFILLABLEGAP) && continue # gap too short to attempt to fill
        lengthsrcgap == lengthtgtgap && lengthsrcgap ≤ MAXIMUMMERGEABLEGAP && continue # gap is going to be merged, don't bother trying to fill
        gapfill!(src2tgt, gap[1], gap[2], alignment, sa, ra, lcps, matchLengthThreshold(lengthsrcgap, lengthtgtgap), mask)
    end

    #= println("links after gapfilling: $(src2tgt.links)")
    for link in src2tgt
        println(link)
    end =#

    # merge adjacent blocks
    link = src2tgt.firstlink
    done = false
        while !done
        if link.next == src2tgt.firstlink; done = true; end
        link = trymergelinks!(src2tgt, link, lenseq1, lenseq2)
    end

    #= println("links after merging: $(src2tgt.links)")
    for link in src2tgt
        println(link)
    end =#

return src2tgt
end

@inline function probMatch(m::Integer, n::Integer, k::Integer)::Float64
    m < k && return 0.0
    n < k && return 0.0
    return 1 - (((1 - 1 / 4^k)^(m - k + 1))^(n - k + 1))
end

function matchLengthThreshold(m::Int32, n::Int32)::Int32
    for k = 1:25
        p = probMatch(m, n, k)
        p < 0.1 && return max(k, 2)
    end
    return 26
end

function gapfill!(mainchain::BlockChain{AlignedBlock}, head::ChainLink{AlignedBlock}, tail::ChainLink{AlignedBlock}, alignment::Alignment, sa::Vector{Int32}, ra::Vector{Int32}, lcps::Vector{Int32}, minblocksize::Int32, mask=nothing)
    #println("gapfilling from $(head.data) to $(tail.data) minblocksize = $(minblocksize)")
    src_length = length(alignment.seq1)
    tgt_length = length(alignment.seq2)
    (srcgap, tgtgap) = contiguousblockgaps(head, tail, src_length, tgt_length)
    chain = blockchain(alignment, sa, ra, lcps, srcgap, tgtgap, minblocksize, mask)
    length(chain) == 0 && return mainchain
    lock(REENTRANT_LOCK)
    head.next = chain.firstlink
    chain.lastlink.next = tail
    mainchain.links += chain.links
    unlock(REENTRANT_LOCK)
    # recursion
    for gap in gaps(head, tail, chain.links + 1)
        (srcgap, tgtgap) = contiguousblockgaps(gap[1], gap[2], src_length, tgt_length)
        lengthsrcgap::Int32 = srcgap.stop - srcgap.start + 1
        lengthtgtgap::Int32 = tgtgap.stop - tgtgap.start + 1
        (lengthsrcgap < MINIMUMFILLABLEGAP || lengthtgtgap < MINIMUMFILLABLEGAP) && continue # gap too short to attempt to fill
        lengthsrcgap == lengthtgtgap && lengthsrcgap ≤ MAXIMUMMERGEABLEGAP && continue # gap is going to be merged, don't bother trying to fill
        mlt = matchLengthThreshold(Int32(length(srcgap)), Int32(length(tgtgap)))
        mlt ≥ minblocksize && continue # gap has already been searched with this minblocksize
        gapfill!(mainchain, gap[1], gap[2], alignment, sa, ra, lcps, mlt, mask)
    end
    return mainchain
end
    
function blockchain(alignment::Alignment, sa::Vector{Int32}, ra::Vector{Int32}, lcps::Vector{Int32}, srcgap, tgtgap, minblocksize, mask=nothing)::BlockChain{AlignedBlock}
    @debug "blockchain from $srcgap to $tgtgap minblocksize = $minblocksize"
    blocks = BlockChain{AlignedBlock}()
    src_length = length(alignment.seq1)
    tgt_length = length(alignment.seq2)
    f = fulcrum(alignment)
    tgt_start = 0
    @inbounds for nt in srcgap
        sa_index = ra[nt]
        offset::Int32 = -1
        minlcp = lcps[sa_index + offset + 1]
        while true
            minlcp < minblocksize && break;
            if sa[sa_index + offset] > f
                tgt_start = sa[sa_index + offset] - f
                tgt_start in tgtgap && break
            end
            offset -= 1
            minlcp = min(minlcp, lcps[sa_index + offset + 1])
        end
        toplcp = minlcp
        toptgt_start = tgt_start
        offset = 1
        tgt_start = 0
        mlt = max(toplcp, minblocksize)
        minlcp = 0
        while true
            if sa_index + offset > length(lcps); minlcp = 0; break; end
            minlcp = minlcp == 0 ? lcps[sa_index + offset] : min(minlcp, lcps[sa_index + offset])
            minlcp < mlt && break;
            if sa[sa_index + offset] > f
                tgt_start = sa[sa_index + offset] - f
                tgt_start in tgtgap && break
            end
            offset += 1
        end
        toplcp = max(toplcp, minlcp)
        toplcp < minblocksize && continue
        if toplcp > minlcp; tgt_start = toptgt_start; end
        new_block = AlignedBlock(mod1(nt, src_length), mod1(tgt_start, tgt_length), mod1(toplcp, min(src_length, tgt_length)))
        if isnothing(mask) || sum(mask[new_block.tgt_index:new_block.tgt_index + new_block.blocklength - 1]) < new_block.blocklength ## add block if it is not included in low complexity region
            append!(blocks, new_block)
        end
    end
    @debug "found $(length(blocks)) blocks"
    return blocks
end

function target_coverage(ff::BlockChain, rf::BlockChain, tgt_length::Integer)
    tgt_coverage = falses(tgt_length)
    for link in ff
        block = link.data
        tgt_coverage[block.tgt_index:min((block.tgt_index + block.blocklength - 1), tgt_length)] .= true
        if block.tgt_index + block.blocklength - 1 > tgt_length
            tgt_coverage[1:mod1(block.tgt_index + block.blocklength - 1, tgt_length)] .= true
        end
    end
    for link in rf
        block = link.data
        tgt_coverage[block.tgt_index:min((block.tgt_index + block.blocklength - 1), tgt_length)] .= true
        if block.tgt_index + block.blocklength - 1 > tgt_length
            tgt_coverage[1:mod1(block.tgt_index + block.blocklength - 1, tgt_length)] .= true
        end
    end
    return sum(tgt_coverage) / tgt_length
end

function entropy_mask(cs::CircularSequence, window_length::Int32)::CircularMask

    comp = zeros(Int32, 4)
    mask = CircularMask(falses(length(cs)))
    entropy = 0.0
    for i::Int32 in 1:window_length
        if cs[i] == DNA_A
            comp[1] += 1
        elseif cs[i] == DNA_C
            comp[2] += 1
        elseif cs[i] == DNA_G
            comp[3] += 1
        elseif cs[i] == DNA_T
            comp[4] += 1
        end
    end
    
    for v in comp
        prob = v / window_length
        prob == 0 && continue
        entropy -= prob * log2(prob)
    end
    entropy ≤ 0 && setindex!(mask, true, 1:window_length)
    
    for i::Int32 in 1:length(cs) - 1
        if cs[i] == DNA_A
            comp[1] -= 1
        elseif cs[i] == DNA_C
            comp[2] -= 1
        elseif cs[i] == DNA_G
            comp[3] -= 1
        elseif cs[i] == DNA_T
            comp[4] -= 1
        end
        if cs[i + window_length] == DNA_A
            comp[1] += 1
        elseif cs[i + window_length] == DNA_C
            comp[2] += 1
        elseif cs[i + window_length] == DNA_G
            comp[3] += 1
        elseif cs[i + window_length] == DNA_T
            comp[4] += 1
        end
        entropy = 0.0
        for v in comp
            prob = v / window_length
            prob == 0 && continue
            entropy -= prob * log2(prob)
        end
        # entropy 0 = homopolymer run
        entropy ≤ 0 && setindex!(mask, true, i:i + window_length - 1)
    end
    
    return mask
end
