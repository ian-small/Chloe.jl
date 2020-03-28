
AlignedBlock = Tuple{Int32,Int32,Int32}
AlignedBlocks = Array{AlignedBlock}

DNAString = AbstractString

function compareSubStrings(a::SubString, b::SubString)
    # a, b = Iterators.Stateful(a), Iterators.Stateful(b)
    count::Int32 = 0
    for (c, d) in zip(a, b)
        c ≠ d && return ifelse(c < d, (count, -1), (count, 1))
        count += 1
    end
    isempty(a) && return ifelse(isempty(b), (count, 0), (count, -1))
    return (count, 1)
end

function alignSAs(a::DNAString, saa::SuffixArray, b::DNAString, sab::SuffixArray)::AlignedBlocks
    lcps = AlignedBlocks(undef, length(saa))
    bpointer = 1
    for apointer = 1:length(saa)
        ssa = SubString(a, saa[apointer])
        oldtuple = (0, 0, 0)
        while true
            ssb = SubString(b, sab[bpointer])
            (lcp, direction) = compareSubStrings(ssa, ssb)
            newtuple = (saa[apointer], sab[bpointer], lcp)
            if (direction < 0 || bpointer == length(sab))
                if (oldtuple[1] == 0) && (bpointer > 1)
                    (lcp2, direction) = compareSubStrings(ssa, SubString(b, sab[bpointer - 1]))
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

function probMatch(m, n, k)
    m < k && return 0
    n < k && return 0
    return 1 - (((1 - 1 / 4^k)^(m - k + 1))^(n - k + 1))
end

function matchLengthThreshold(m, n)
    for k = 1:25
        p = probMatch(m, n, k)
        p < 0.1 && return max(k, 2)
    end
    return 26
end

function lcps2AlignmentBlocks(lcps::AlignedBlocks, circular::Bool, min_run_length::Integer)::AlignedBlocks
    seqlen = length(lcps)
    aligned_blocks = AlignedBlocks(undef, 0)
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
    sum(map(b->b[3], blocks))
end

function fillGap(block1::AlignedBlock, block2::AlignedBlock, 
        refloop::DNAString, refSA::SuffixArray, refRA::SuffixArray, 
        targetloop::DNAString, targetSA::SuffixArray, targetRA::SuffixArray)::AlignedBlocks

    ref_gap = block2[1] - block1[1] - block1[3]
    target_gap = block2[2] - block1[2] - block1[3]

    ref_gap < 5 && return []
    target_gap < 5 && return []
    ref_gap_range = block1[1] + block1[3]:block2[1] - 1
    target_gap_range = block1[2] + block1[3]:block2[2] - 1

    # make gap SAs from genome SAs
    ref_ranks_slice = sort(refRA[ref_gap_range])
    target_ranks_slice = sort(targetRA[target_gap_range])

    ref_gap_SA = SuffixArray(undef, ref_gap)
    target_gap_SA = SuffixArray(undef, target_gap)
    for i = 1:ref_gap
        ref_gap_SA[i] = refSA[ref_ranks_slice[i]]
    end
    for i = 1:target_gap
        target_gap_SA[i] = targetSA[target_ranks_slice[i]]
    end

    # align gap SAs to get lcps
    gap_lcps = alignSAs(refloop, ref_gap_SA, targetloop, target_gap_SA)
    gap_blocks = lcps2AlignmentBlocks(gap_lcps, false, matchLengthThreshold(ref_gap, target_gap))
    return gap_blocks
end
MaybeAlignedBlock = Union{AlignedBlock,Nothing}
function mergeBlocks(block1::AlignedBlock, block2::AlignedBlock)::Tuple{AlignedBlock,MaybeAlignedBlock}
    gap1 = block2[1] - block1[1] - block1[3]
    gap2 = block2[2] - block1[2] - block1[3]
    if gap1 == gap2 && gap1 < 12
        p = probMatch(block1[3] + gap1, block2[3] + gap2, block2[3])
        # println(block1," ",block2," ",p)
        if p < 0.01
            return (block1[1], block1[2], block1[3] + gap1 + block2[3]), nothing
        end
    end
    return block1, block2
end

function fillAllGaps!(aligned_blocks::AlignedBlocks, 
    refloop::DNAString, refSA::SuffixArray, refRA, 
    targetloop::DNAString, targetSA::SuffixArray, targetRA::SuffixArray)::AlignedBlocks
    block1_pointer = 1
    block2_pointer = 2
    while block2_pointer <= length(aligned_blocks)
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
        new_blocks = fillGap(aligned_blocks[block1_pointer], aligned_blocks[block2_pointer], 
                refloop, refSA, refRA, targetloop, targetSA, targetRA)
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
    for i = 1:length(blocks)
        block = blocks[i]
        rev_blocks[i] = (a_length - block[3] - block[1] + 2, b_length - block[3] - block[2] + 2, block[3])
    end
    return rev_blocks
end

function mergeBlockArrays(blocks1::AlignedBlocks, blocks2::AlignedBlocks)::AlignedBlocks
    # assume sorted arrays
    merged_array = AlignedBlocks(undef, 0)
    blocks1_pointer = 1
    blocks2_pointer = 1
    while blocks1_pointer <= length(blocks1) && blocks2_pointer <= length(blocks2)
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
    while blocks1_pointer <= length(blocks1)
        push!(merged_array, blocks1[blocks1_pointer])
        blocks1_pointer += 1
    end
    while blocks2_pointer <= length(blocks2)
        block2 = blocks2[blocks2_pointer]
        push!(merged_array, (block2[2], block2[1], block2[3]))
        blocks2_pointer += 1
    end
    return merged_array
end

AlignedBlock = Tuple{Int32,Int32,Int32}
AlignedBlocks = Array{AlignedBlock}

function alignLoops(refloop::DNAString, 
                    refSA::SuffixArray, refRA::SuffixArray, 
                    targetloop::DNAString,
                    target_SA::SuffixArray, target_RA::SuffixArray)::Tuple{AlignedBlocks,AlignedBlocks}
    # check if sequences are identical, allowing for rotation
    if length(refloop) == length(targetloop)
        match = findfirst(SubString(refloop, 1, floor(Int, (length(refloop) + 1) / 2)), targetloop)
        if !isnothing(match)
            aligned_block = [(1, match[1], length(refloop))]
            return aligned_block, revCompBlocks(aligned_block, length(refSA), length(target_SA))
        end
    end

    # align from ref to target
    lcps = alignSAs(refloop, refSA, targetloop, target_SA)
    # set minimum match length to 18 based on matchLengthThreshold() calculation for 50-150kbp sequences
    rt_aligned_blocks = lcps2AlignmentBlocks(lcps, true, matchLengthThreshold(length(refSA), length(target_SA)))
    rt_aligned_blocks = fillAllGaps!(rt_aligned_blocks, refloop, refSA, refRA, targetloop, target_SA, target_RA)
    # print("Coverage ref to target: ")
    # println(blockCoverage(rt_aligned_blocks))

    # align from target to ref
    lcps = alignSAs(targetloop, target_SA, refloop, refSA)
    # set minimum match length to 18 based on matchLengthThreshold() calculation for 50-150kbp sequences
    tr_aligned_blocks = lcps2AlignmentBlocks(lcps, true, matchLengthThreshold(length(refSA), length(target_SA)))
    tr_aligned_blocks = fillAllGaps!(tr_aligned_blocks, targetloop, target_SA, target_RA, refloop, refSA, refRA)
    # print("Coverage target to ref: ")
    # println(blockCoverage(tr_aligned_blocks))


    # combine rt and tr blocks; tr blocks need [1] and [2] to be swapped
    # for tr_block in tr_aligned_blocks
    #     push!(rt_aligned_blocks,(tr_block[2],tr_block[1],tr_block[3]))
    # end
    merged_blocks = mergeBlockArrays(rt_aligned_blocks, tr_aligned_blocks)

    # print("Merged coverage: ")
    # println(blockCoverage(merged_blocks))

    rev_blocks = revCompBlocks(merged_blocks, length(refSA), length(target_SA))
    # println("Coverage: ", blockCoverage(merged_blocks) + blockCoverage(rev_blocks))
    return merged_blocks, rev_blocks
end
