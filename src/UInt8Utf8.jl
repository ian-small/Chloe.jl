import Base

ncodeunits(s::Vector{UInt8}) = length(s)
codeunit(s::Vector{UInt8}) = UInt8

codeunit(s::Vector{UInt8}, i::Int) = s[i]

# Ugh! I had to copy and paste this from julia/base/strings/string.jl
# why didn't they put the UTF8 logic on a Vector{UInt8}? (which
# is basically what a String is)

## thisind, nextind ##

Base.@propagate_inbounds thisind(s::Vector{UInt8}, i::Int) = thisind_mm(s, i)


@inline function thisind_mm(s::Vector{UInt8}, i::Int)
    i == 0 && return 0
    n = ncodeunits(s)
    i == n + 1 && return i
    @boundscheck Base.between(i, 1, n) || throw(Base.BoundsError(s, i))
    @inbounds b = codeunit(s, i)
    (b & 0xc0 == 0x80) & (i - 1 > 0) || return i
    @inbounds b = codeunit(s, i - 1)
    Base.between(b, 0b11000000, 0b11110111) && return i - 1
    (b & 0xc0 == 0x80) & (i - 2 > 0) || return i
    @inbounds b = codeunit(s, i - 2)
    Base.between(b, 0b11100000, 0b11110111) && return i - 2
    (b & 0xc0 == 0x80) & (i - 3 > 0) || return i
    @inbounds b = codeunit(s, i - 3)
    Base.between(b, 0b11110000, 0b11110111) && return i - 3
    return i
end

Base.@propagate_inbounds nextind(s::Vector{UInt8}, i::Int) = nextind_mm(s, i)


@inline function nextind_mm(s::Vector{UInt8}, i::Int)
    i == 0 && return 1
    n = ncodeunits(s)
    @boundscheck Base.between(i, 1, n) || throw(Base.BoundsError(s, i))
    @inbounds l = codeunit(s, i)
    (l < 0x80) | (0xf8 ≤ l) && return i + 1
    if l < 0xc0
        i′ = thisind(s, i)
        return i′ < i ? nextind(s, i′) : i + 1
    end
    # first continuation byte
    (i += 1) > n && return i
    @inbounds b = codeunit(s, i)
    b & 0xc0 ≠ 0x80 && return i
    ((i += 1) > n) | (l < 0xe0) && return i
    # second continuation byte
    @inbounds b = codeunit(s, i)
    b & 0xc0 ≠ 0x80 && return i
    ((i += 1) > n) | (l < 0xf0) && return i
    # third continuation byte
    @inbounds b = codeunit(s, i)
    ifelse(b & 0xc0 ≠ 0x80, i, i + 1)
end



isvalid(s::Vector{UInt8}) = Base.isvalid(String, s)

## required core functionality ##

Base.@propagate_inbounds function iterate(s::Vector{UInt8}, i::Int=firstindex(s))
    i > ncodeunits(s) && return nothing
    b = codeunit(s, i)
    u = UInt32(b) << 24
    Base.between(b, 0x80, 0xf7) || return reinterpret(Char, u), i + 1
    return iterate_continued_mm(s, i, u)
end



function iterate_continued_mm(s::Vector{UInt8}, i::Int, u::UInt32)
    u < 0xc0000000 && (i += 1; @goto ret)
    n = ncodeunits(s)
    # first continuation byte
    (i += 1) > n && @goto ret
    @inbounds b = codeunit(s, i)
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b) << 16
    # second continuation byte
    ((i += 1) > n) | (u < 0xe0000000) && @goto ret
    @inbounds b = codeunit(s, i)
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b) << 8
    # third continuation byte
    ((i += 1) > n) | (u < 0xf0000000) && @goto ret
    @inbounds b = codeunit(s, i)
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b); i += 1
    @label ret
    return reinterpret(Char, u), i
end

Base.@propagate_inbounds function getindex(s::Vector{UInt8}, i::Int)
    b = codeunit(s, i)
    u = UInt32(b) << 24
    Base.between(b, 0x80, 0xf7) || return reinterpret(Char, u)
    return getindex_continued(s, i, u)
end


function getindex_continued(s::Vector{UInt8}, i::Int, u::UInt32)
    if u < 0xc0000000
        # called from `getindex` which checks bounds
        @inbounds isvalid(s, i) && @goto ret
        Base.string_index_err(s, i)
    end
    n = ncodeunits(s)

    (i += 1) > n && @goto ret
    @inbounds b = codeunit(s, i) # cont byte 1
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b) << 16

    ((i += 1) > n) | (u < 0xe0000000) && @goto ret
    @inbounds b = codeunit(s, i) # cont byte 2
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b) << 8

    ((i += 1) > n) | (u < 0xf0000000) && @goto ret
    @inbounds b = codeunit(s, i) # cont byte 3
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b)
    @label ret
    return reinterpret(Char, u)
end

getindex(s::Vector{UInt8}, r::UnitRange{<:Integer}) = s[Int(first(r)):Int(last(r))]

@inline function getindex(s::Vector{UInt8}, r::UnitRange{Int})
    isempty(r) && return ""
    i, j = first(r), last(r)
    @boundscheck begin
        checkbounds(s, r)
        @inbounds isvalid(s, i) || Base.string_index_err(s, i)
        @inbounds isvalid(s, j) || Base.string_index_err(s, j)
    end
    j = nextind(s, j) - 1
    n = j - i + 1
    ss = Base._string_n(n)
    GC.@preserve s ss Base.unsafe_copyto!(pointer(ss), pointer(s, i), n)
    return ss
end

utf8_length(s::Vector{UInt8}) = length_continued(s, 1, ncodeunits(s), ncodeunits(s))

@inline function utf8_length(s::Vector{UInt8}, i::Int, j::Int)
    @boundscheck begin
        0 < i ≤ ncodeunits(s) + 1 || throw(Base.BoundsError(s, i))
        0 ≤ j < ncodeunits(s) + 1 || throw(Base.BoundsError(s, j))
    end
    j < i && return 0
    @inbounds i, k = thisind(s, i), i
    c = j - i + (i == k)
    length_continued(s, i, j, c)
end


@inline function length_continued(s::Vector{UInt8}, i::Int, n::Int, c::Int)
    i < n || return c
    @inbounds b = codeunit(s, i)
    @inbounds while true
        while true
            (i += 1) ≤ n || return c
            0xc0 ≤ b ≤ 0xf7 && break
            b = codeunit(s, i)
        end
        l = b
        b = codeunit(s, i) # cont byte 1
        c -= (x = b & 0xc0 == 0x80)
        x & (l ≥ 0xe0) || continue

        (i += 1) ≤ n || return c
        b = codeunit(s, i) # cont byte 2
        c -= (x = b & 0xc0 == 0x80)
        x & (l ≥ 0xf0) || continue

        (i += 1) ≤ n || return c
        b = codeunit(s, i) # cont byte 3
        c -= (b & 0xc0 == 0x80)
    end
end

## overload methods for efficiency ##

isvalid(s::Vector{UInt8}, i::Int) = checkbounds(Bool, s, i) && thisind(s, i) == i

function isascii(s::Vector{UInt8})
    @inbounds for i = 1:ncodeunits(s)
        codeunit(s, i) >= 0x80 && return false
    end
    return true
end

