
AlignedBlock = Tuple{Int32,Int32,Int32}
AlignedBlocks = Vector{AlignedBlock}

import .MappedString: MMappedString, ASCII

datasize(a::AlignedBlocks) = length(a) * sizeof(AlignedBlock)

@inline function compareSubStrings(a::SubString{MMappedString{ASCII}}, 
                                   b::SubString{MMappedString{ASCII}})::Tuple{Int32,Int}
    la, lb = length(a), length(b)
    m = min(la, lb)

    count::Int32 = 0
    one::Int32 = 1
    ia = a.string.ptr
    ib = b.string.ptr
    # compare bytes!
    @inbounds for i in 1:m
        c = ia[i + a.offset]
        d = ib[i + b.offset]
        c ≠ d && return c < d ? (count, -1) : (count, 1)
        count += one
    end
    la == 0 && return lb == 0 ? (count, 0) : (count, -1)
    return (count, 1)
end

@inline function compareSubStrings(a::SubString, b::SubString)::Tuple{Int32,Int}
    # a, b = Iterators.Stateful(a), Iterators.Stateful(b)
    count::Int32 = 0
    one::Int32 = 1
    for (c, d) in zip(a, b)
        c ≠ d && return c < d ? (count, -1) : (count, 1)
        count += one
    end
    isempty(a) && return isempty(b) ? (count, 0) : (count, -1)
    return (count, 1)
end

function alignSAs(src::DNAString, srcSA::SuffixArray, tgt::DNAString, tgtSA::SuffixArray)::AlignedBlocks
    # src is the reference suffix position
    # tgt is the suffixes of what we want to align
    lcps = AlignedBlocks(undef, length(srcSA))
    tgt_idx = 1
    zerot = (zero(Int32), zero(Int32), zero(Int32))
    @inbounds for src_idx = 1:length(srcSA)
        ssa = SubString(src, srcSA[src_idx])
        oldtuple = zerot
        while true
            ssb = SubString(tgt, tgtSA[tgt_idx])
            (lcp, direction) = compareSubStrings(ssa, ssb)
            newtuple = (srcSA[src_idx], tgtSA[tgt_idx], lcp)
            if (direction < 0 || tgt_idx == length(tgtSA))
                if (oldtuple[1] == 0) && (tgt_idx > 1)
                    sab1 = SubString(tgt, tgtSA[tgt_idx - 1])
                    (lcp2, direction) = compareSubStrings(ssa, sab1)
                    oldtuple = (srcSA[src_idx], tgtSA[tgt_idx - 1], lcp2)
                end
                if (lcp >= oldtuple[3])
                    lcps[src_idx] = newtuple
                else
                    lcps[src_idx] = oldtuple
                end
                break
            end
            tgt_idx += 1
            oldtuple = newtuple
        end
    end
    return sort!(lcps)
end

@inline function probMatch(m::Integer, n::Integer, k::Integer)::Float64
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

function old_lcps2AlignmentBlocks(lcps::AlignedBlocks, circular::Bool, min_run_length::Integer)::AlignedBlocks
    # misses check of last element of lcps
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
function lcps2AlignmentBlocks(lcps::AlignedBlocks, circular::Bool, min_run_length::Integer)::AlignedBlocks
    # length of aligned_blocks is usually much small than lcps
    aligned_blocks = AlignedBlocks()

    b_start = -1

    for lcp in lcps
        if lcp[2] == b_start + 1 # target location
            b_start += 1
            continue
        end

        if lcp[3] >= min_run_length
            if isempty(aligned_blocks)
                push!(aligned_blocks, lcp)
            else
                # see if we overlap with last added block; if so, only keep longest block
                p_start, _, p_run_length = aligned_blocks[end]
                a_start, _, a_run_length = lcp
                if rangesOverlap(a_start, a_run_length, p_start, p_run_length)
                    if a_run_length > p_run_length
                        aligned_blocks[end] = lcp
                    end
                else
                    push!(aligned_blocks, lcp)
                end
            end
        end
        b_start = lcp[2]

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

    lcps2AlignmentBlocks(gap_lcps, false, matchLengthThreshold(ref_gap, tgt_gap))
end


MaybeAlignedBlock = Union{AlignedBlock,Nothing}
@inline function mergeBlocks(block1::AlignedBlock, block2::AlignedBlock)::Tuple{AlignedBlock,MaybeAlignedBlock}
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
    tgtloop::DNAString, tgtSA::SuffixArray, tgtRA::SuffixArray)
    block1_index = 1
    block2_index = 2
    
    # nshifted should be **small**
    # @debug
    start_len = length(aligned_blocks)
    nsplices = nshifted = ndeleted = nnew = 0
    
    @inbounds while block2_index <= length(aligned_blocks)
        block1 = aligned_blocks[block1_index]
        block2 = aligned_blocks[block2_index]
        block1, block2 = mergeBlocks(block1, block2)
        while isnothing(block2) # blocks1 and 2 got merged

            # @debug
            nsplices += 1
            nshifted = length(aligned_blocks) - block2_index
            ndeleted = + block2_index - block1_index

            splice!(aligned_blocks, block1_index:block2_index, [block1]) # replace block1 and block2 with merged block
            if block2_index > length(aligned_blocks)
                return aligned_blocks
            end
            block2 = aligned_blocks[block2_index]
            block1, block2 = mergeBlocks(block1, block2)
        end
        new_blocks = fillGap(aligned_blocks[block1_index],
                             aligned_blocks[block2_index], 
                             refloop, refSA, refRA,
                             tgtloop, tgtSA, tgtRA)
        if (isempty(new_blocks))
            block1_index += 1
            block2_index += 1
        else
            # @debug
            nsplices += 1
            nshifted = length(aligned_blocks) - block2_index + 1
            nnew += length(new_blocks)

            # insert before block index 2
            splice!(aligned_blocks, block2_index:block2_index - 1, new_blocks)
        end
    end
    @debug "fillAllGaps: splices=$(nsplices)  gaps=$(nnew) shifted=$(nshifted) deleted=$(ndeleted) len=$(start_len) -> $(length(aligned_blocks))"
end

function revCompBlocks(blocks::AlignedBlocks, a_length::Integer, b_length::Integer)::AlignedBlocks
    rev_blocks = AlignedBlocks(undef, length(blocks))
    two, a_length, b_length = Int32(2), Int32(a_length), Int32(b_length)
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

    blocks1_index = 1
    blocks2_index = 1
    blocks1_len = length(blocks1)
    blocks2_len = length(blocks2)

    merged_array = AlignedBlocks()
    # this is basically true for the tests I've used....
    sizehint!(merged_array, blocks1_len + blocks2_len)

    @inbounds while blocks1_index <= blocks1_len && blocks2_index <= blocks2_len
        block1 = blocks1[blocks1_index]
        block2 = blocks2[blocks2_index]
        if block1[1] == block2[2] && block1[2] == block2[1] && block1[3] == block2[3] # both the same, only add one of them
            push!(merged_array, block1)
            blocks1_index += 1
            blocks2_index += 1
        elseif block1[1] + block1[3] < block2[2] # block1 can't overlap anything in blocks2, so add it
            push!(merged_array, block1)
            blocks1_index += 1
        elseif block2[2] + block2[3] < block1[1] # block2 can't overlap anything in blocks1, so add it
            push!(merged_array, (block2[2], block2[1], block2[3]))
            blocks2_index += 1
        # blocks must overlap if we got this far
        elseif block1[1] - block1[2] == block2[2] - block2[1] # same offset so compatible
            blocka_start = min(block1[1], block2[2])
            blockb_start = min(block1[2], block2[1])
            block_length = max(block1[1] + block1[3], block2[2] + block2[3]) - blocka_start + one(Int32)
            push!(merged_array, (blocka_start, blockb_start, block_length)) # add merger of the two blocks
            blocks1_index += 1
            blocks2_index += 1
        else # overlap but incompatible, add both
            push!(merged_array, block1)
            push!(merged_array, (block2[2], block2[1], block2[3]))
            blocks1_index += 1
            blocks2_index += 1
        end
    end
    @inbounds while blocks1_index <= blocks1_len
        push!(merged_array, blocks1[blocks1_index])
        blocks1_index += 1
    end
    @inbounds while blocks2_index <= blocks2_len
        block2 = blocks2[blocks2_index]
        push!(merged_array, (block2[2], block2[1], block2[3]))
        blocks2_index += 1
    end
    @debug "mergeBlockArrays: sizehint! $(length(merged_array)) ≅ $(blocks1_len + blocks2_len)"
    return merged_array
end

function alignHelper(src_id::String, strand::Char,
    src::DNAString, srcSA::SuffixArray, srcRA::SuffixArray,
    tgt::DNAString, tgtSA::SuffixArray, tgtRA::SuffixArray)
    # ~30% for alignSAs and 60% for fillAllGaps!
    lcps = alignSAs(src, srcSA, tgt, tgtSA)
    aligned_blocks = lcps2AlignmentBlocks(lcps, true, matchLengthThreshold(length(srcSA), length(tgtSA)))
    fillAllGaps!(aligned_blocks, src, srcSA, srcRA, tgt, tgtSA, tgtRA)
    @debug "alignLoops: [$(src_id)]$(strand) lcps#=$(length(lcps)) -> aligned#=$(length(aligned_blocks))"
    aligned_blocks
end

function alignLoops(src_id::String,
                    ref_loop::DNAString, ref_SA::SuffixArray, ref_RA::SuffixArray, 
                    tgt_loop::DNAString, tgt_SA::SuffixArray, tgt_RA::SuffixArray;
                    rev=true)::Tuple{AlignedBlocks,Union{Nothing,AlignedBlocks}}
    
    # check if sequences are identical, allowing for rotation
    nr, nt = length(ref_loop), length(tgt_loop)
    if nr == nt
        match = findfirst(SubString(ref_loop, 1, floor(Int, (nr + 1) / 2)), tgt_loop)
        if !isnothing(match)
            aligned_blocks = [(one(Int32), Int32(first(match)), Int32(nr))]
            return aligned_blocks, rev ? revCompBlocks(aligned_blocks, length(ref_SA), length(tgt_SA)) : nothing
        end
    end

    function rt()
        alignHelper(src_id, '+', ref_loop, ref_SA, ref_RA, tgt_loop, tgt_SA, tgt_RA)
    end
    function tr()
        alignHelper(src_id, '-', tgt_loop, tgt_SA, tgt_RA, ref_loop, ref_SA, ref_RA)
    end

    blockf, blockr = fetch.((Threads.@spawn w()) for w in [rt, tr])

    merged_blocks = mergeBlockArrays(blockf, blockr)
    return merged_blocks, rev ? revCompBlocks(merged_blocks, length(ref_SA), length(tgt_SA)) : nothing
end
