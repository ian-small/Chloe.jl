using FASTX
include("sais.jl")
include("ChainedBlocks.jl")

global const MINIMUMFILLABLEGAP = 20
global const MAXIMUMMERGEABLEGAP = 400
global const one = Int32(1)

function align(; fasta_files=Vector{String}, directory=nothing)
    reader = open(FASTA.Reader, fasta_files[1])
    record = FASTA.Record()
    read!(reader, record)
    seq1 = FASTA.sequence(record)
    close(reader)
    reader = open(FASTA.Reader, fasta_files[2])
    record = FASTA.Record()
    read!(reader, record)
    seq2 = FASTA.sequence(record)
    close(reader)
    
    local mainchain
    for i in 1:11
        @time align2strings(seq1,seq2)
    end

    # open("At_vs_Nt.sais.merged.txt", "w") do outfile
    #     # println(length(mainchain))
    #     for link in mainchain
    #         write(outfile,string(link),'\n')
    #     end
    # end
    #coverage(mainchain, length(seq1), length(seq2))
end

@inline function is_source(index::Int32, fulcrum::Int32)
    return index < fulcrum ? true : false
end

@inline function is_target(index::Int32, fulcrum::Int32)
    return index > fulcrum ? true : false
end

function align2seqs(seq1::LongDNASeq,seq2::LongDNASeq)::BlockChain{AlignedBlock}
    #extend and join strings; currently assumes both sequence are circular
    #TO DO: generalise for either or both of seq1, seq2 as linear sequences
    alignment = Alignment(seq1, seq2)
    sa = sais(alignment, zeros(Int32,length(alignment)), 0, length(alignment), 16) .+= 1 #shift to 1-based julia indexes
    ra = invperm(sa)
    lcps = lcparray(alignment, sa, ra)
    lenseq1 = Int32(length(seq1))
    lenseq2 = Int32(length(seq2))

    # sa_index = ra[1]
    # for i in -5:5
    #     println(sa[sa_index+i], '\t', lcps[sa_index+i], '\t', joinedstring[sa[sa_index+i]:min(sa[sa_index+i]+50,end)])
    # end

    src2tgt = blockchain(alignment, sa, ra, lcps, one:lenseq1, one:lenseq2, matchLengthThreshold(lenseq1, lenseq2))
    src2tgt = circularise(src2tgt, lenseq1)

    @debug "links before gapfilling: $(src2tgt.links)"
    @debug begin
        for link in src2tgt
            println(link)
        end
    end

    Threads.@threads for gap in gaps(src2tgt)
        (srcgap, tgtgap) = contiguousblockgaps(gap[1], gap[2], lenseq1, lenseq2)
        length(srcgap) < MINIMUMFILLABLEGAP && length(tgtgap) < MINIMUMFILLABLEGAP && continue #gap too short to attempt to fill
        length(srcgap) == length(tgtgap) && length(srcgap) ≤ MAXIMUMMERGEABLEGAP && continue #gap is going to be merged, don't bother trying to fill
        gapfill!(src2tgt, gap[1], gap[2], alignment, sa, ra, lcps, matchLengthThreshold(Int32(length(srcgap)), Int32(length(tgtgap))))
    end

    @debug "links after gapfilling: $(src2tgt.links)"
    @debug begin
        for link in src2tgt
            println(link)
        end
    end

    #merge adjacent blocks
    link = src2tgt.firstlink
    done = false
    while !done
        if link.next == src2tgt.firstlink; done = true; end
        link = trymergelinks!(src2tgt, link, lenseq1, lenseq2)
    end

    @debug "links after merging: $(src2tgt.links)"
    @debug begin
        for link in src2tgt
            println(link)
        end
    end

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

function gapfill!(mainchain::BlockChain{AlignedBlock}, head::ChainLink{AlignedBlock}, tail::ChainLink{AlignedBlock}, alignment::Alignment, sa::Vector{Int32}, ra::Vector{Int32}, lcps::Vector{Int32}, minblocksize::Int32)
    @debug "gapfilling from $(head.data) to $(tail.data) minblocksize = $(minblocksize)"
    src_length = length(alignment.seq1)
    tgt_length = length(alignment.seq2)
    (srcgap, tgtgap) = contiguousblockgaps(head, tail, src_length, tgt_length)
    chain = blockchain(alignment, sa, ra, lcps, srcgap, tgtgap, minblocksize)
    # println("chain length: ", length(chain))
    # for link in chain
    #     println(link)
    # end
    length(chain) == 0 && return mainchain
    lock(REENTRANT_LOCK)
    head.next = chain.firstlink
    chain.lastlink.next = tail
    mainchain.links += chain.links
    unlock(REENTRANT_LOCK)
    #recursion
    for gap in gaps(head,tail,chain.links+1)
        #println("maybe fill gap: ", gap)
        (srcgap, tgtgap) = contiguousblockgaps(gap[1], gap[2], src_length, tgt_length)
        length(srcgap) < MINIMUMFILLABLEGAP || length(tgtgap) < MINIMUMFILLABLEGAP && continue #gap too short to attempt to fill
        length(srcgap) == length(tgtgap) && length(srcgap) ≤ MAXIMUMMERGEABLEGAP && continue #gap is going to be merged, don't bother trying to fill
        mlt = matchLengthThreshold(Int32(length(srcgap)), Int32(length(tgtgap)))
        #println("mlt: ", mlt)
        mlt ≥ minblocksize && continue #gap has already been searched with this minblocksize
        gapfill!(mainchain, gap[1], gap[2], alignment, sa, ra, lcps, mlt)
    end
    return mainchain
end
    
function blockchain(alignment::Alignment, sa::Vector{Int32}, ra::Vector{Int32}, lcps::Vector{Int32}, srcgap, tgtgap, minblocksize)::BlockChain{AlignedBlock}
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
            minlcp = min(minlcp,lcps[sa_index + offset + 1])
        end
        toplcp = minlcp
        toptgt_start = tgt_start
        #println(nt,"\toffset: ", offset, "\ttoplowlcp: ", toplcp)
        offset = 1
        tgt_start = 0
        mlt = max(toplcp, minblocksize)
        minlcp = 0
        while true
            if sa_index + offset > length(lcps); minlcp = 0; break; end
            minlcp = minlcp == 0 ? lcps[sa_index + offset] : min(minlcp,lcps[sa_index + offset])
            minlcp < mlt && break;
            if sa[sa_index + offset] > f
                tgt_start = sa[sa_index + offset] - f
                tgt_start in tgtgap && break
            end
            offset += 1
            #println("offset: ",offset,"\tmincp: ", minlcp)
        end
        toplcp = max(toplcp,minlcp)
        toplcp < minblocksize && continue
        if toplcp > minlcp; tgt_start = toptgt_start; end
        new_block = AlignedBlock(mod1(nt, src_length), mod1(tgt_start, tgt_length), mod1(toplcp, min(src_length, tgt_length)))
        append!(blocks, new_block)
        # if new_block.blocklength > 10000
        #     @error "excessive blocklength in blockchain() $(new_block)"
        #     println(nt,"\toffset: ", offset, "\tsa_index: ", sa_index, "\tsa length: ", length(sa), "\tsa: ", sa[sa_index + offset], "\ttgt_start: ", tgt_start, "\tmlt: ", mlt, "\ttoplcp: ", toplcp)
        #     for i in -5:5
        #         println(sa[sa_index+i], '\t', lcps[sa_index+i], '\t', LongDNASeq(reinterpret(DNA, alignment[sa[sa_index+i]:min(sa[sa_index+i]+50,end)])) )
        #     end
        # end
    end
    @debug "found $(length(blocks)) blocks"
    return blocks
end

function target_coverage(ff, rf, tgt_length::Integer)
    tgt_coverage = falses(tgt_length)
    for block in ff
        tgt_coverage[block.tgt_index:min((block.tgt_index + block.blocklength - 1),tgt_length)] .= true
        if block.tgt_index + block.blocklength - 1 > tgt_length
            tgt_coverage[1:mod1(block.tgt_index + block.blocklength - 1, tgt_length)] .= true
        end
    end
    for block in rf
        tgt_coverage[block.tgt_index:min((block.tgt_index + block.blocklength - 1),tgt_length)] .= true
        if block.tgt_index + block.blocklength - 1 > tgt_length
            tgt_coverage[1:mod1(block.tgt_index + block.blocklength - 1, tgt_length)] .= true
        end
    end
    return sum(tgt_coverage)/tgt_length
end
