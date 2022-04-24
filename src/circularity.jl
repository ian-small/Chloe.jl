using BioSequences

const KMERSIZE = 10
const ENTROPYTHRESHOLD = 0
const SKETCHSIZE = 400

struct CircularVector
    v::Vector{Int32}
end

@inline Base.length(cv::CircularVector) = Int32(length(cv.v))
@inline Base.getindex(cv::CircularVector, i::Int32) = @inbounds getindex(cv.v, mod1(i, length(cv)))
@inline Base.getindex(cv::CircularVector, is::Vector{<:Integer}) = getindex(cv.v, mod1.(is,length(cv.v)))
function Base.getindex(cv::CircularVector, r::UnitRange{<:Integer})
    len = length(cv)
    start = mod1(r.start, len)
    stop = mod1(r.stop, len)
    if start < stop
        return cv.v[start:stop]
    else
        return vcat(cv.v[start:end], cv.v[1:stop])
    end
end
@inline Base.setindex!(cv::CircularVector, value::Int32, i::Int32) = @inbounds setindex!(cv.v, value, mod1(i, length(cv)))
@inline Base.push!(cv::CircularVector, value::Int32) = @inbounds push!(cv.v, value)
@inline Base.sum(cv::CircularVector) = sum(cv.v)
function Base.sum(cv::CircularVector, r::UnitRange{<:Integer})
    sum = 0
    for i in r
        sum += cv[i]
    end
    return sum
end

function circularintersect(r1::UnitRange, r2::UnitRange, d::Int32)::UnitRange{Int32}
    maxintersect = UnitRange{Int32}(0, -1)
    maxintersectlength = 0
    i1 = intersect(r1, r2)
    if length(i1) > maxintersectlength
        maxintersect = i1
        maxintersectlength = length(i1)
    end
    if r2.stop > d
        i2 = intersect(r1 .+ d, r2)
        if length(i2) > maxintersectlength
            maxintersect = i2
            maxintersectlength = length(i2)
        end
    end
    if r1.stop > d
        i3 = intersect(r1, r2 .+ d)
        if length(i3) > maxintersectlength
            maxintersect = i3
            maxintersectlength = length(i3)
        end
    end
    if maxintersect.start > d
        maxintersect = maxintersect .- d
    end
    return maxintersect
end

struct CircularMask
    m::BitVector
end

@inline Base.length(cm::CircularMask) = Int32(length(cm.m))
@inline Base.getindex(cm::CircularMask, i::Int32) = @inbounds getindex(cm.m, mod1(i, length(cm)))
function Base.getindex(cm::CircularMask, r::UnitRange{<:Integer})
    len = length(cm)
    start = mod1(r.start, len)
    stop = mod1(r.stop, len)
    if start < stop
        return cm.m[start:stop]
    else
        return vcat(cm.m[start:end], cm.m[1:stop])
    end
end
@inline Base.setindex!(cm::CircularMask, value::Bool, i::Int32) = @inbounds setindex!(cm.m, value, mod1(i, length(cm)))
function Base.setindex!(cm::CircularMask, value::Bool, r::UnitRange{<:Integer})
    for i in r
        setindex!(cm.m, value, mod1(i, length(cm)))
    end
end
function Base.sum(cm::CircularMask, r::UnitRange{<:Integer})
    len = length(cm)
    start = mod1(r.start, len)
    stop = mod1(r.stop, len)
    if start < stop
        return sum(cm.m[start:stop])
    else
        return sum(cm.m[start:end]) + sum(cm.m[1:stop])
    end
end
Base.iterate(cm::CircularMask) = iterate(cm.m)
Base.iterate(cm::CircularMask, state) = iterate(cm.m, state)
Base.reverse(cm::CircularMask) = CircularMask(reverse(cm.m))

function entropy_mask!(seq::LongDNA{4}, mask::CircularMask)
    comp, entropy = sequence_entropy(seq[1:KMERSIZE])
    entropy ≤ 0 && setindex!(mask, true, 1:KMERSIZE)
    for i::Int32 in 1:length(seq)-KMERSIZE
        if seq[i] == DNA_A
            comp[1] -= 1
        elseif seq[i] == DNA_C
            comp[2] -= 1
        elseif seq[i] == DNA_G
            comp[3] -= 1
        elseif seq[i] == DNA_T
            comp[4] -= 1
        end
        if seq[i+KMERSIZE] == DNA_A
            comp[1] += 1
        elseif seq[i+KMERSIZE] == DNA_C
            comp[2] += 1
        elseif seq[i+KMERSIZE] == DNA_G
            comp[3] += 1
        elseif seq[i+KMERSIZE] == DNA_T
            comp[4] += 1
        end
        entropy = 0.0
        for v in comp
            prob = v / KMERSIZE
            prob == 0 && continue
            entropy -= prob * log2(prob)
        end
        # entropy 0 = homopolymer run
        entropy ≤ 0 && setindex!(mask, true, i:i+KMERSIZE-1)
    end
end

struct CircularSequence
    length::Int32
    sequence::LongDNA{2}
    mask::CircularMask  #used to mask low entropy regions, ambiguous bases, gaps

    function CircularSequence(seq::LongDNA{2}, mask::CircularMask)
        new(length(seq), append!(seq, LongSubSeq(seq, 1:length(seq)-1)), mask)
    end

    function CircularSequence(seq::LongDNA{4})
        mask = CircularMask(falses(length(seq)))
        entropy_mask!(seq, mask)
        for (n::Int32, nt) in enumerate(seq)
            if isambiguous(nt) || isgap(nt)
                setindex!(mask, true, n)
                seq[n] = DNA_A
            end
        end
        compressedseq = LongDNA{2}(seq)
        new(length(compressedseq), append!(compressedseq, LongSubSeq(compressedseq, 1:length(compressedseq)-1)), mask)
    end
end

@inline Base.length(cs::CircularSequence) = cs.length
@inline Base.getindex(cs::CircularSequence, i::Int32) = @inbounds getindex(cs.sequence, mod1(i, cs.length))

function Base.getindex(cs::CircularSequence, r::UnitRange{<:Integer})
    @assert length(r) <= length(cs)
    if r.start > length(cs) || r.start < 1
        r = range(mod1(r.start, cs.length); length=length(r))
    end
    return LongSubSeq(cs.sequence, r)
end

function reverse_complement(cs::CircularSequence)
    rc = BioSequences.reverse_complement(cs.sequence[1:cs.length])
    return CircularSequence(rc, reverse(cs.mask))
end

@inline function getcodon(cs::CircularSequence, index::Int32)
    return (cs[index], cs[index+Int32(1)], cs[index+Int32(2)])
end