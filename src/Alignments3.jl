
AlignedBlock = Tuple{Int32,Int32,Int32}
# AlignedBlock = @NamedTuple begin
#     src::Int32
#     tgt::Int32
#     length::Int32
# end
# import Base
# Base.convert(::Type{AlignedBlock}, t::Tuple{Int32,Int32,Int32}) = AlignedBlock((t[1], t[2], t[3]))
# Base.convert(::Type{AlignedBlock}, t::Tuple{Int32,Int32,Int64}) = AlignedBlock((t[1], t[2], t[3]))

# @inline src_end(b::AlignedBlock) = b.src + b.length
# @inline tgt_end(b::AlignedBlock) = b.tgt + b.length

const Two = Int32(2)
const Three = Int32(3)


AlignedBlocks = Vector{AlignedBlock}


@inline function compareSubStrings(a::SubString, b::SubString)::Tuple{Int32,Int}
    # a, b = Iterators.Stateful(a), Iterators.Stateful(b)
    count::Int32 = 0
    for (c, d) in zip(a, b)
        c ≠ d && return c < d ? (count, -1) : (count, 1)
        count += 1
    end
    isempty(a) && return isempty(b) ? (count, 0) : (count, -1)
    return (count, 1)
end

function alignSAs(a::DNAString, saa::SuffixArray, b::DNAString, sab::SuffixArray)::AlignedBlocks
    lcps = AlignedBlocks(undef, length(saa))
    bpointer = 1
    zerot = (zero(Int32), zero(Int32), zero(Int32))
    @inbounds for apointer = 1:length(saa)
        ssa = SubString(a, saa[apointer])
        oldtuple = zerot
        while true
            ssb = SubString(b, sab[bpointer])
            (lcp, direction) = compareSubStrings(ssa, ssb)
            newtuple = (saa[apointer], sab[bpointer], lcp)
            if (direction < 0 || bpointer == length(sab))
                if (oldtuple[1] == 0) && (bpointer > 1)
                    sab1 = SubString(b, sab[bpointer - 1])
                    (lcp2, direction) = compareSubStrings(ssa, sab1)
                    oldtuple = (saa[apointer], sab[bpointer - 1], lcp2)
                end
                if (lcp >= oldtuple[3])
                    lcps[apointer] = newtuple
                else
                    lcps[apointer] = oldtuple
                end
                break
            end
            bpointer += 1
            oldtuple = newtuple
        end
    end
    return sort!(lcps)
end

function probMatch(m::Integer, n::Integer, k::Integer)::Float64
    m < k && return 0.0
    n < k && return 0.0
    return 1 - (((1 - 1 / 4^k)^(m - k + 1))^(n - k + 1))
end

function matchLengthThreshold(m::T, n::T)::Int where {T <: Integer}
    for k = 1:25
        p = probMatch(m, n, k)
        p < 0.1 && return max(k, 2)
    end
    return 26
end

function lcps2AlignmentBlocks(lcps::AlignedBlocks, circular::Bool, min_run_length::Integer)::AlignedBlocks
    seqlen = length(lcps)
    aligned_blocks = AlignedBlocks()
    a_start = lcps[1][1]
    b_start = lcps[1][2]
    pointer = -1
    run_length = 0
    for lcp in lcps
        if lcp[2] == pointer + 1
            pointer += 1
        else
            if run_length >= min_run_length
                # see if we overlap with last added block; if so, only keep longest block
                if !isempty(aligned_blocks) && rangesOverlap(a_start, run_length, last(aligned_blocks)[1], last(aligned_blocks)[3])
                    if run_length > last(aligned_blocks)[3]
                        pop!(aligned_blocks)
                        push!(aligned_blocks, (a_start, b_start, run_length))
                    end
                else
                    push!(aligned_blocks, (a_start, b_start, run_length))
                end
            end
            run_length = lcp[3]
            a_start = lcp[1]
            b_start = lcp[2]
            pointer = b_start
        end
    end
    return aligned_blocks
end

function blockCoverage(blocks::AlignedBlocks)
    sum(map(b -> b[3], blocks))
end


function fillGap(block1::AlignedBlock, block2::AlignedBlock, 
        refloop::DNAString, refSA::SuffixArray, refRA::SuffixArray, 
        tgtloop::DNAString, tgtSA::SuffixArray, tgtRA::SuffixArray)::AlignedBlocks

    ref_gap = block2[1] - block1[1] - block1[3]
    tgt_gap = block2[2] - block1[2] - block1[3]

    ref_gap < 5 && return []
    tgt_gap < 5 && return []
    ref_gap_range = block1[1] + block1[3]:block2[1] - 1
    tgt_gap_range = block1[2] + block1[3]:block2[2] - 1

    # make gap SAs from genome SAs
    # REM: v[l:h] makes a *copy* so we can do an inplace sort (of the copy)
    ref_ranks_slice = sort!(refRA[ref_gap_range])
    tgt_ranks_slice = sort!(tgtRA[tgt_gap_range])

    ref_gap_SA = SuffixArray(undef, ref_gap)
    tgt_gap_SA = SuffixArray(undef, tgt_gap)
    @inbounds for i = 1:ref_gap
        ref_gap_SA[i] = refSA[ref_ranks_slice[i]]
    end
    @inbounds for i = 1:tgt_gap
        tgt_gap_SA[i] = tgtSA[tgt_ranks_slice[i]]
    end
 
    # align gap SAs to get lcps
    gap_lcps = alignSAs(refloop, ref_gap_SA, tgtloop, tgt_gap_SA)
    gap_blocks = lcps2AlignmentBlocks(gap_lcps, false, matchLengthThreshold(ref_gap, tgt_gap))
    return gap_blocks
end


MaybeAlignedBlock = Union{AlignedBlock,Nothing}
function mergeBlocks(block1::AlignedBlock, block2::AlignedBlock)::Tuple{AlignedBlock,MaybeAlignedBlock}
    gap1 = block2[1] - block1[1] - block1[3]
    gap2 = block2[2] - block1[2] - block1[3]
    if gap1 == gap2 && gap1 < 12
        p = probMatch(block1[3] + gap1, block2[3] + gap2, block2[3])
        if p < 0.01
            return (block1[1], block1[2], block1[3] + gap1 + block2[3]), nothing
        end
    end
    return block1, block2
end

# this is the killer....
# most of Chloë is spent in fillGap
function fillAllGaps!(aligned_blocks::AlignedBlocks, 
    refloop::DNAString, refSA::SuffixArray, refRA::SuffixArray, 
    tgtloop::DNAString, tgtSA::SuffixArray, tgtRA::SuffixArray)::AlignedBlocks
    block1_pointer = 1
    block2_pointer = 2
    @inbounds while block2_pointer <= length(aligned_blocks)
        block1 = aligned_blocks[block1_pointer]
        block2 = aligned_blocks[block2_pointer]
        block1, block2 = mergeBlocks(block1, block2)
        while isnothing(block2) # blocks1 and 2 got merged
            splice!(aligned_blocks, block1_pointer:block2_pointer, [block1]) # replace block1 and block2 with merged block
            if block2_pointer > length(aligned_blocks)
                return aligned_blocks
            end
            block2 = aligned_blocks[block2_pointer]
            block1, block2 = mergeBlocks(block1, block2)
        end
        new_blocks = fillGap(aligned_blocks[block1_pointer],
                             aligned_blocks[block2_pointer], 
                             refloop, refSA, refRA,
                             tgtloop, tgtSA, tgtRA)
        if (isempty(new_blocks))
            block1_pointer += 1
            block2_pointer += 1
        else
            splice!(aligned_blocks, block2_pointer:block2_pointer - 1, new_blocks)
        end
    end
    return aligned_blocks
end

function revCompBlocks(blocks::AlignedBlocks, a_length::Integer, b_length::Integer)::AlignedBlocks
    rev_blocks = AlignedBlocks(undef, length(blocks))
    two, a_length, b_length = Two, Int32(a_length), Int32(b_length)
    @inbounds for i = 1:length(blocks)
        block = blocks[i]
        rev_blocks[i] = (a_length - block[3] - block[1] + two, 
                         b_length - block[3] - block[2] + two,
                         block[3])
    end
    return rev_blocks
end

function mergeBlockArrays(blocks1::AlignedBlocks, blocks2::AlignedBlocks)::AlignedBlocks
    # assume sorted arrays
    merged_array = AlignedBlocks()

    blocks1_pointer = 1
    blocks2_pointer = 1
    blocks1_len = length(blocks1)
    blocks2_len = length(blocks2)

    sizehint!(merged_array, blocks1_len + blocks2_len)

    @inbounds while blocks1_pointer <= blocks1_len && blocks2_pointer <= blocks2_len
        block1 = blocks1[blocks1_pointer]
        block2 = blocks2[blocks2_pointer]
        if block1[1] == block2[2] && block1[2] == block2[1] && block1[3] == block2[3] # both the same, only add one of them
            push!(merged_array, block1)
            blocks1_pointer += 1
            blocks2_pointer += 1
        elseif block1[1] + block1[3] < block2[2] # block1 can't overlap anything in blocks2, so add it
            push!(merged_array, block1)
            blocks1_pointer += 1
        elseif block2[2] + block2[3] < block1[1] # block2 can't overlap anything in blocks1, so add it
            push!(merged_array, (block2[2], block2[1], block2[3]))
            blocks2_pointer += 1
        # blocks must overlap if we got this far
        elseif block1[1] - block1[2] == block2[2] - block2[1] # same offset so compatible
            blocka_start = min(block1[1], block2[2])
            blockb_start = min(block1[2], block2[1])
            block_length = max(block1[1] + block1[3], block2[2] + block2[3]) - blocka_start + 1
            push!(merged_array, (blocka_start, blockb_start, block_length)) # add merger of the two blocks
            blocks1_pointer += 1
            blocks2_pointer += 1
        else # overlap but incompatible, add both
            push!(merged_array, block1)
            push!(merged_array, (block2[2], block2[1], block2[3]))
            blocks1_pointer += 1
            blocks2_pointer += 1
        end
    end
    @inbounds while blocks1_pointer <= blocks1_len
        push!(merged_array, blocks1[blocks1_pointer])
        blocks1_pointer += 1
    end
    @inbounds while blocks2_pointer <= blocks2_len
        block2 = blocks2[blocks2_pointer]
        push!(merged_array, (block2[2], block2[1], block2[3]))
        blocks2_pointer += 1
    end
    @debug "mergeBlockArrays: $(length(merged_array)) ≅ $(blocks1_len + blocks2_len)"
    return merged_array
end

function alignLoops(ref_loop::DNAString, ref_SA::SuffixArray, ref_RA::SuffixArray, 
                    tgt_loop::DNAString, tgt_SA::SuffixArray, tgt_RA::SuffixArray)::Tuple{AlignedBlocks,AlignedBlocks}
    # check if sequences are identical, allowing for rotation
    if length(ref_loop) == length(tgt_loop)
        match = findfirst(SubString(ref_loop, 1, floor(Int, (length(ref_loop) + 1) / 2)), tgt_loop)
        if !isnothing(match)
            aligned_blocks = [(one(Int32), Int32(match[1]), Int32(length(ref_loop)))]
            return aligned_blocks, revCompBlocks(aligned_blocks, length(ref_SA), length(tgt_SA))
        end
    end
    function align(src::DNAString, srcSA::SuffixArray, srcRA::SuffixArray,
                   tgt::DNAString, tgtSA::SuffixArray, tgtRA::SuffixArray)
        # ~30% for alignSAs and 60% for fillAllGaps!
        lcps = alignSAs(src, srcSA, tgt, tgtSA)
        @debug "align[$(Threads.threadid())]" lcps = length(lcps)
        aligned_blocks = lcps2AlignmentBlocks(lcps, true, matchLengthThreshold(length(srcSA), length(tgtSA)))
        aligned_blocks = fillAllGaps!(aligned_blocks, src, srcSA, srcRA, tgt, tgtSA, tgtRA)
        aligned_blocks
    end
    block = Vector{AlignedBlocks}(undef, 2)
    function rt()
        block[1] = align(ref_loop, ref_SA, ref_RA, tgt_loop, tgt_SA, tgt_RA)
    end
    function tr()
        block[2] = align(tgt_loop, tgt_SA, tgt_RA, ref_loop, ref_SA, ref_RA)
    end

    Threads.@threads for worker in [rt, tr]
        worker()
    end

    merged_blocks = mergeBlockArrays(block[1], block[2])

    rev_blocks = revCompBlocks(merged_blocks, length(ref_SA), length(tgt_SA))
    return merged_blocks, rev_blocks
end
