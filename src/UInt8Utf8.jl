
module UInt8Utf8

export utf8_codeunit, utf8_ncodeunits, utf8_getindex, utf8_isascii, utf8_isvalid
export utf8_iterate, utf8_length, utf8_nextind, utf8_prevind, utf8_thisind

export ascii_getindex, ascii_isvalid, ascii_length, ascii_iterate
export ascii_nextind, ascii_prevind, ascii_thisind

import Base

utf8_ncodeunits(s::Vector{UInt8}) = length(s)
utf8_codeunit(s::Vector{UInt8}) = UInt8
@inline utf8_codeunit(s::Vector{UInt8}, i::Int) = s[i]

# Ugh! I had to copy and paste this from julia/base/strings/string.jl
# why didn't they put the UTF8 logic on a Vector{UInt8}? (which
# is basically what a String is)

## thisind, nextind ##

Base.@propagate_inbounds utf8_thisind(s::Vector{UInt8}, i::Int) = thisind_mm(s, i)

@inline function thisind_mm(s::Vector{UInt8}, i::Int)
    i == 0 && return 0
    n = utf8_ncodeunits(s)
    i == n + 1 && return i
    @boundscheck Base.between(i, 1, n) || throw(Base.BoundsError(s, i))
    @inbounds b = utf8_codeunit(s, i)
    (b & 0xc0 == 0x80) & (i - 1 > 0) || return i
    @inbounds b = utf8_codeunit(s, i - 1)
    Base.between(b, 0b11000000, 0b11110111) && return i - 1
    (b & 0xc0 == 0x80) & (i - 2 > 0) || return i
    @inbounds b = utf8_codeunit(s, i - 2)
    Base.between(b, 0b11100000, 0b11110111) && return i - 2
    (b & 0xc0 == 0x80) & (i - 3 > 0) || return i
    @inbounds b = utf8_codeunit(s, i - 3)
    Base.between(b, 0b11110000, 0b11110111) && return i - 3
    return i
end

Base.@propagate_inbounds utf8_nextind(s::Vector{UInt8}, i::Int) = nextind_mm(s, i)

@inline function nextind_mm(s::Vector{UInt8}, i::Int)
    i == 0 && return 1
    n = utf8_ncodeunits(s)
    @boundscheck Base.between(i, 1, n) || throw(Base.BoundsError(s, i))
    @inbounds l = utf8_codeunit(s, i)
    (l < 0x80) | (0xf8 ≤ l) && return i + 1
    if l < 0xc0
        i′ = utf8_thisind(s, i)
        return i′ < i ? utf8_nextind(s, i′) : i + 1
    end
    # first continuation byte
    (i += 1) > n && return i
    @inbounds b = utf8_codeunit(s, i)
    b & 0xc0 ≠ 0x80 && return i
    ((i += 1) > n) | (l < 0xe0) && return i
    # second continuation byte
    @inbounds b = utf8_codeunit(s, i)
    b & 0xc0 ≠ 0x80 && return i
    ((i += 1) > n) | (l < 0xf0) && return i
    # third continuation byte
    @inbounds b = utf8_codeunit(s, i)
    ifelse(b & 0xc0 ≠ 0x80, i, i + 1)
end

utf8_prevind(s::Vector{UInt8}, i::Int) = utf8_prevind(s, i, 1)

function utf8_prevind(s::Vector{UInt8}, i::Int, n::Int)
    n < 0 && throw(Base.ArgumentError("n cannot be negative: $n"))
    z = utf8_ncodeunits(s) + 1
    @boundscheck 0 < i ≤ z || throw(BoundsError(s, i))
    n == 0 && return utf8_thisind(s, i) == i ? i : Base.string_index_err(s, i)
    while n > 0 && 1 < i
        @inbounds n -= utf8_isvalid(s, i -= 1)
    end
    return i - n
end

utf8_isvalid(s::Vector{UInt8}) = Base.isvalid(String, s)

## required core functionality ##

Base.@propagate_inbounds function utf8_iterate(s::Vector{UInt8}, i::Int)
    i > utf8_ncodeunits(s) && return nothing
    b = utf8_codeunit(s, i)
    u = UInt32(b) << 24
    Base.between(b, 0x80, 0xf7) || return reinterpret(Char, u), i + 1
    return iterate_continued_mm(s, i, u)
end

function iterate_continued_mm(s::Vector{UInt8}, i::Int, u::UInt32)
    u < 0xc0000000 && (i += 1; @goto ret)
    n = utf8_ncodeunits(s)
    # first continuation byte
    (i += 1) > n && @goto ret
    @inbounds b = utf8_codeunit(s, i)
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b) << 16
    # second continuation byte
    ((i += 1) > n) | (u < 0xe0000000) && @goto ret
    @inbounds b = utf8_codeunit(s, i)
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b) << 8
    # third continuation byte
    ((i += 1) > n) | (u < 0xf0000000) && @goto ret
    @inbounds b = utf8_codeunit(s, i)
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b); i += 1
    @label ret
    return reinterpret(Char, u), i
end

Base.@propagate_inbounds function utf8_getindex(s::Vector{UInt8}, i::Int)
    b = utf8_codeunit(s, i)
    u = UInt32(b) << 24
    Base.between(b, 0x80, 0xf7) || return reinterpret(Char, u)
    return getindex_continued(s, i, u)
end

function getindex_continued(s::Vector{UInt8}, i::Int, u::UInt32)
    if u < 0xc0000000
        # called from `getindex` which checks bounds
        @inbounds utf8_isvalid(s, i) && @goto ret
        Base.string_index_err(s, i)
    end
    n = utf8_ncodeunits(s)

    (i += 1) > n && @goto ret
    @inbounds b = utf8_codeunit(s, i) # cont byte 1
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b) << 16

    ((i += 1) > n) | (u < 0xe0000000) && @goto ret
    @inbounds b = utf8_codeunit(s, i) # cont byte 2
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b) << 8

    ((i += 1) > n) | (u < 0xf0000000) && @goto ret
    @inbounds b = utf8_codeunit(s, i) # cont byte 3
    b & 0xc0 == 0x80 || @goto ret
    u |= UInt32(b)
    @label ret
    return reinterpret(Char, u)
end

@inline function utf8_getindex(s::Vector{UInt8}, r::UnitRange{Int})
    isempty(r) && return ""
    i, j = first(r), last(r)
    @boundscheck begin
        checkbounds(s, r)
        @inbounds utf8_isvalid(s, i) || Base.string_index_err(s, i)
        @inbounds utf8_isvalid(s, j) || Base.string_index_err(s, j)
    end
    j = utf8_nextind(s, j) - 1
    n = j - i + 1
    ss = Base._string_n(n)
    GC.@preserve s ss Base.unsafe_copyto!(pointer(ss), pointer(s, i), n)
    return ss
end

utf8_length(s::Vector{UInt8}) = length_continued(s, 1, utf8_ncodeunits(s), utf8_ncodeunits(s))

@inline function utf8_length(s::Vector{UInt8}, i::Int, j::Int)
    @boundscheck begin
        0 < i ≤ utf8_ncodeunits(s) + 1 || throw(Base.BoundsError(s, i))
        0 ≤ j < utf8_ncodeunits(s) + 1 || throw(Base.BoundsError(s, j))
    end
    j < i && return 0
    @inbounds i, k = utf8_thisind(s, i), i
    c = j - i + (i == k)
    length_continued(s, i, j, c)
end

@inline function length_continued(s::Vector{UInt8}, i::Int, n::Int, c::Int)
    i < n || return c
    @inbounds b = utf8_codeunit(s, i)
    @inbounds while true
        while true
            (i += 1) ≤ n || return c
            0xc0 ≤ b ≤ 0xf7 && break
            b = utf8_codeunit(s, i)
        end
        l = b
        b = utf8_codeunit(s, i) # cont byte 1
        c -= (x = b & 0xc0 == 0x80)
        x & (l ≥ 0xe0) || continue

        (i += 1) ≤ n || return c
        b = utf8_codeunit(s, i) # cont byte 2
        c -= (x = b & 0xc0 == 0x80)
        x & (l ≥ 0xf0) || continue

        (i += 1) ≤ n || return c
        b = utf8_codeunit(s, i) # cont byte 3
        c -= (b & 0xc0 == 0x80)
    end
end

## overload methods for efficiency ##

utf8_isvalid(s::Vector{UInt8}, i::Int) = checkbounds(Bool, s, i) && utf8_thisind(s, i) == i

function utf8_isascii(s::Vector{UInt8})
    @inbounds for i = 1:utf8_ncodeunits(s)
        utf8_codeunit(s, i) >= 0x80 && return false
    end
    return true
end

# assuming ascii

Base.@propagate_inbounds function ascii_iterate(s::Vector{UInt8}, i::Int)
    i > utf8_ncodeunits(s) && return nothing
    b = utf8_codeunit(s, i)
    Char(b), i + 1
end

# step forward and back by one byte == one character
Base.@propagate_inbounds function ascii_nextind(s::Vector{UInt8}, i::Int)
    i == 0 && return 1
    n = utf8_ncodeunits(s)
    @boundscheck Base.between(i, 1, n) || throw(Base.BoundsError(s, i))
    i + 1
end

Base.@propagate_inbounds function ascii_thisind(s::Vector{UInt8}, i::Int)
    i == 0 && return 0
    n = utf8_ncodeunits(s)
    i == n + 1 && return i
    @boundscheck Base.between(i, 1, n) || throw(Base.BoundsError(s, i))
    i
end

Base.@propagate_inbounds function ascii_prevind(s::Vector{UInt8}, i::Int)
    i == 1 && return 0
    n = utf8_ncodeunits(s)
    @boundscheck Base.between(i, 1, n + 1) || throw(Base.BoundsError(s, i))
    i - 1
end

@inline function ascii_length(s::Vector{UInt8}, i::Int, j::Int)
    @boundscheck begin
        m = utf8_ncodeunits(s)
        0 < i ≤ m + 1 || throw(Base.BoundsError(s, i))
        0 ≤ j < m + 1 || throw(Base.BoundsError(s, j))
    end
    j < i && return 0
    return j - i + 1
end

ascii_isvalid(s::Vector{UInt8}, i::Int) = checkbounds(Bool, s, i)

Base.@propagate_inbounds function ascii_getindex(s::Vector{UInt8}, i::Int)
    b = utf8_codeunit(s, i)
    Char(b)
end

Base.@propagate_inbounds function ascii_getindex(s::Vector{UInt8}, r::UnitRange{Int})
    isempty(r) && return ""
    return (@view s[r]) .|> Char |> String
end

end

