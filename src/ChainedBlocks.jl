struct AlignedBlock <: AbstractInterval{Int32}
    src_index::Int32
    tgt_index::Int32
    blocklength::Int32
end

Base.first(b::AlignedBlock) = b.src_index
Base.last(b::AlignedBlock) = b.src_index + b.blocklength - one(Int32)

const emptyblock = AlignedBlock(zero(Int32), zero(Int32), zero(Int32))

mutable struct ChainLink{AlignedBlock}
    data::AlignedBlock
    next::ChainLink{AlignedBlock}
    function ChainLink{AlignedBlock}()
        link = new(emptyblock)
        link.next = link
        return link
    end
    function ChainLink{AlignedBlock}(data::AlignedBlock)
        link = new(data)
        link.next = link
        return link
    end
end

mutable struct BlockChain{AlignedBlock}
    links::Integer
    firstlink::ChainLink{AlignedBlock}
    lastlink::ChainLink{AlignedBlock}
    function BlockChain{AlignedBlock}()
        chain = new()
        link = ChainLink{AlignedBlock}()
        chain.links = 0
        chain.firstlink = link
        chain.lastlink = link
        return chain
    end
    function BlockChain{AlignedBlock}(data::AlignedBlock)
        chain = new()
        link = ChainLink{AlignedBlock}(data)
        chain.links = 1
        chain.firstlink = link
        chain.lastlink = link
        return chain
    end
end

function Base.length(chain::BlockChain{AlignedBlock})
    return chain.links
end

function Base.isempty(chain::BlockChain{AlignedBlock})
    return chain.links == 0 ? true : false
end

function Base.isempty(link::ChainLink{AlignedBlock})
    return link.data === emptyblock ? true : false
end

function Base.string(link::ChainLink{AlignedBlock})
    return join([string(link.data.src_index),string(link.data.tgt_index),string(link.data.blocklength)], " ")
end

function Base.show(io::IO, link::ChainLink{AlignedBlock})
    return print(link)
end

function Base.print(io::IO, link::ChainLink{AlignedBlock})
    return print(string(link))
end

function Base.append!(chain::BlockChain{AlignedBlock}, b::AlignedBlock)
    # check if the blocks overlap before appending
    if chain.lastlink.data.src_index + chain.lastlink.data.blocklength ≥ b.src_index
        if chain.lastlink.data.blocklength ≥ b.blocklength
            return chain # new block is worse than existing link, so don't append
        else
            if isempty(chain); chain.links = 1; end
            chain.lastlink.data = b # new block is better than existing link, so replace it
        end
    else # append new link
        link = ChainLink{AlignedBlock}(b)
        if isempty(chain)
            chain.firstlink = link
        end
        chain.lastlink.next = link
        chain.lastlink = link
        chain.links += 1
    end
    return chain
end

function Base.append!(mainchain::BlockChain{AlignedBlock}, appendage::BlockChain{AlignedBlock})
    if isempty(chain)
        mainchain.firstlink = appendage.firstlink
    end
    mainchain.lastlink.next = appendage.firstlink
    mainchain.lastlink = appendage.lastlink
    mainchain.links += appendage.links
    return mainchain
end

function contiguousblockgaps(link1::ChainLink{AlignedBlock}, link2::ChainLink{AlignedBlock}, src_length, tgt_length)

    srcgapstart::Int32 = 1
    tgtgapstart::Int32 = 1
    srcgapend::Int32 = src_length
    tgtgapend::Int32 = tgt_length

    if !isempty(link1)
        block1 = link1.data
        srcgapstart = mod1(block1.src_index + block1.blocklength, src_length)
        tgtgapstart = mod1(block1.tgt_index + block1.blocklength, tgt_length)
    end
    if !isempty(link2)
        block2 = link2.data
        srcgapend = mod1(block2.src_index - 1, src_length)
        tgtgapend = mod1(block2.tgt_index - 1, tgt_length)
    end
    return (srcgapstart:srcgapend, tgtgapstart:tgtgapend)
end

function circularise(chain::BlockChain{AlignedBlock}, genome_length::Int32)
    # last link may overlap first link; if so remove firstlink and link to firstlink.next
    if mod1(chain.lastlink.data.src_index + chain.lastlink.data.blocklength, genome_length) == chain.firstlink.data.src_index + chain.firstlink.data.blocklength
        chain.lastlink.next = chain.firstlink.next
        chain.firstlink = chain.lastlink.next
        chain.links -= 1
    else
        chain.lastlink.next = chain.firstlink
    end
    return chain
end

function trymergelinks!(chain::BlockChain{AlignedBlock}, link1::ChainLink{AlignedBlock}, lenseq1, lenseq2)
    link2 = link1.next
    link1 == link2 && return link1
    (srcgap, tgtgap) = contiguousblockgaps(link1, link2, lenseq1, lenseq2)
    lengthsrcgap::Int32 = srcgap.stop - srcgap.start + 1
    lengthtgtgap::Int32 = tgtgap.stop - tgtgap.start + 1
    if lengthsrcgap ==  lengthtgtgap && 0 < lengthsrcgap ≤ MAXIMUMMERGEABLEGAP
        #println("merging: ",link1,"\t",link2, "\t", srcgap, "\t", tgtgap)
        # safe to not mod1() the indexes as they don't change
        link1.data = AlignedBlock(link1.data.src_index, link1.data.tgt_index, link1.data.blocklength + length(srcgap) + link2.data.blocklength)
        link1.next = link2.next
        if link2 == chain.lastlink
            chain.lastlink = link1
        end
        if link2 == chain.firstlink
            chain.firstlink = link1.next
        end
        chain.links -= 1
        return link1
    end
    return link2
end

# Iteration

@inline function Base.iterate(chain::BlockChain{AlignedBlock})
    return (chain.firstlink, chain.firstlink)
end

@inline function Base.iterate(chain::BlockChain{AlignedBlock}, link::ChainLink{AlignedBlock})
    link.next == link && return nothing # end of linear chain
    link.next == chain.firstlink && return nothing # end of circular chain
    return (link.next, link.next)
end

function gaps(chain::BlockChain{AlignedBlock})::Vector{Tuple{ChainLink{AlignedBlock},ChainLink{AlignedBlock}}}
    circular = chain.lastlink.next == chain.firstlink ? true : false
    if circular; gapslength = length(chain)
    else; gapslength = length(chain) + 1; end
    gaps = Vector{Tuple{ChainLink{AlignedBlock},ChainLink{AlignedBlock}}}(undef, gapslength)
    gapslength == 0 && return gaps
    # set internal gaps
    link1 = chain.firstlink
    for i in 2:chain.links
        link2 = link1.next
        gaps[i] = (link1, link2)
        link1 = link2
    end
    # set terminal gaps
    if circular
        gaps[1] = (chain.lastlink, chain.firstlink);
    else
        gaps[1] = (ChainLink{AlignedBlock}(), chain.firstlink);
        gaps[chain.links + 1] = (chain.lastlink, ChainLink{AlignedBlock}())
    end
    return gaps
end

function gaps(head::ChainLink{AlignedBlock}, tail::ChainLink{AlignedBlock}, numgaps)::Vector{Tuple{ChainLink{AlignedBlock},ChainLink{AlignedBlock}}}
    gaps = Vector{Tuple{ChainLink{AlignedBlock},ChainLink{AlignedBlock}}}(undef, numgaps)
    link1 = head
    for i in 1:numgaps
        link2 = link1.next
        gaps[i] = (link1, link2)
        link1 = link2
    end
    @assert link1 == tail
    return gaps
end

AlignedBlocks = Vector{AlignedBlock}
datasize(a::AlignedBlocks) = length(a) * sizeof(AlignedBlock)

function ll2vector(chain::BlockChain{AlignedBlock})::AlignedBlocks
    v = AlignedBlocks(undef,chain.links)
    length(v) == 0 && return v
    for (i,link) in enumerate(chain)
        v[i] = link.data
    end
    return v
end

struct BlockTree
    btree::IntervalTree{Int32,AlignedBlock}
    length::Int32
    wrapped_intervals::Dict{AlignedBlock,AlignedBlock}
end

function Base.push!(ctree::BlockTree, interval::AlignedBlock)
    if last(interval) ≤ ctree.length
        push!(ctree.btree, interval)
    else
        first_interval = AlignedBlock(first(interval), zero(Int32), ctree.length - first(interval) + one(Int32))
        second_interval = AlignedBlock(1, zero(Int32), last(interval) - ctree.length)
        push!(ctree.btree, first_interval)
        ctree.wrapped_intervals[first_interval] = interval
        push!(ctree.btree, second_interval)
        ctree.wrapped_intervals[second_interval] = interval
    end
end

function chain2tree(chain::BlockChain{AlignedBlock}, glength::Int32)::BlockTree
    blocktree = BlockTree(IntervalTree{Int32,AlignedBlock}(), glength, Dict{Interval{Int32},AlignedBlock}())
    for link in chain
        push!(blocktree, link.data)
    end
    return blocktree
end

function chain2rctree(chain::BlockChain{AlignedBlock}, src_length::Int32, tgt_length::Int32)::BlockTree
    blocktree = BlockTree(IntervalTree{Int32,AlignedBlock}(), src_length, Dict{Interval{Int32},AlignedBlock}())
    for link in chain
        b = link.data
        push!(blocktree, AlignedBlock(mod1(src_length - b.blocklength - b.src_index + Int32(2), src_length), mod1(tgt_length - b.blocklength - b.tgt_index + Int32(2), tgt_length), b.blocklength))
    end
    return blocktree
end