
## Abstract String backed by a readonly Memory Mapped Vector{UInt8}
# as I can't be sure that String(Vector{UInt8}) doesn't make a copy.... :(
struct Unicode end
struct ASCII end

struct MMappedString{T <: Union{Unicode,ASCII}} <: AbstractString
    ptr::Vector{UInt8}
end

MMappedString(ptr) = MMappedString{Unicode}(ptr)

MMappedString(s::String) = MMappedString(Vector{UInt8}(s)) # copy unfortunately
MMappedString{ASCII}(s::String) = MMappedString{ASCII}(Vector{UInt8}(s)) # copy unfortunately
import Base
import Base: ==

function ==(x::MMappedString, y::MMappedString)
    x.ptr == y.ptr
end
function ==(x::String, y::MMappedString)
    sa = sizeof(x)
    sa == length(y.ptr) && Base._memcmp(x, pointer(y.ptr), sa) == 0
end

==(y::MMappedString, x::String) = x == y

@inline function Base.cmp(a::MMappedString, b::MMappedString)
    Base.cmp(a.ptr, b.ptr)
end


function Base.cmp(a::String, b::MMappedString)
    al, bl = sizeof(a), length(b)
    c = Base._memcmp(a, pointer(b.ptr), min(al, bl))
    return c < 0 ? -1 : c > 0 ? +1 : Base.cmp(al, bl)
end

@inline function Base.cmp(a::MMappedString, b::String)
    -Base.cmp(b, a)
end


@inline function Base.hash(s::MMappedString, h::UInt)
    Base.hash(s.ptr, h)
end

function Base.show(io::IO, s::MMappedString)
    print(io, "MMappedString[$(length(s.ptr))] @ $(pointer(s.ptr)) of type $(codeunit(s))")
end

# Base.print(io::IO, s::MMappedString) = print(io, string(s))
# Base.textwidth(s::MMappedString) = textwidth(string(s))


# Base.string(x::MMappedString) = String(x)

# Base.convert(::Type{MMappedString}, x::String) = MMappedString(x)
# Base.convert(::Type{String}, x::MMappedString) = String(x)
# Base.String(x::MMappedString) = String(copy(x.ptr))

Base.pointer(s::MMappedString) = pointer(s.ptr)
Base.pointer(s::MMappedString, i::Integer) = pointer(s.ptr, i)

Base.ncodeunits(s::MMappedString) = length(s.ptr)
Base.codeunit(s::MMappedString) = UInt8

Base.@propagate_inbounds function Base.codeunit(s::MMappedString, i::Int)
    s.ptr[i]
end

# Ugh! I had to copy and paste this from julia/base/strings/string.jl
# why didn't they put the UTF8 logic on a Vector{UInt8}

## thisind, nextind ##

Base.@propagate_inbounds Base.thisind(s::MMappedString, i::Int) = _thisind_str(s, i)


@inline function _thisind_str(s::MMappedString, i::Int)
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

Base.@propagate_inbounds Base.nextind(s::MMappedString, i::Int) = _nextind_str(s, i)


@inline function _nextind_str(s::MMappedString, i::Int)
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



Base.isvalid(s::MMappedString) = Base.isvalid(String, s.ptr)

## required core functionality ##

Base.@propagate_inbounds function Base.iterate(s::MMappedString, i::Int=firstindex(s))
    i > ncodeunits(s) && return nothing
    b = codeunit(s, i)
    u = UInt32(b) << 24
    Base.between(b, 0x80, 0xf7) || return reinterpret(Char, u), i + 1
    return iterate_continued(s, i, u)
end



function iterate_continued(s::MMappedString, i::Int, u::UInt32)
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

Base.@propagate_inbounds function Base.getindex(s::MMappedString, i::Int)
    b = codeunit(s, i)
    u = UInt32(b) << 24
    Base.between(b, 0x80, 0xf7) || return reinterpret(Char, u)
    return getindex_continued(s, i, u)
end


function getindex_continued(s::MMappedString, i::Int, u::UInt32)
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

Base.getindex(s::MMappedString, r::UnitRange{<:Integer}) = s[Int(first(r)):Int(last(r))]

@inline function Base.getindex(s::MMappedString, r::UnitRange{Int})
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

Base.length(s::MMappedString) = length_continued(s, 1, ncodeunits(s), ncodeunits(s))

@inline function Base.length(s::MMappedString, i::Int, j::Int)
    @boundscheck begin
        0 < i ≤ ncodeunits(s) + 1 || throw(Base.BoundsError(s, i))
        0 ≤ j < ncodeunits(s) + 1 || throw(Base.BoundsError(s, j))
    end
    j < i && return 0
    @inbounds i, k = thisind(s, i), i
    c = j - i + (i == k)
    length_continued(s, i, j, c)
end


@inline function length_continued(s::MMappedString, i::Int, n::Int, c::Int)
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

Base.isvalid(s::MMappedString, i::Int) = checkbounds(Bool, s, i) && thisind(s, i) == i

function Base.isascii(s::MMappedString)
    @inbounds for i = 1:ncodeunits(s)
        codeunit(s, i) >= 0x80 && return false
    end
    return true
end


# specializations

## Must be ASCII! Otherwise the bitwise equality above
# would not work and latin1 from 0x80 up are 2 byte representations

function Base.show(io::IO, s::MMappedString{ASCII})
    print(io, "MMappedString[$(length(s.ptr))] @ $(pointer(s.ptr)) of type $(codeunit(s)) Latin1")
end
Base.@propagate_inbounds function Base.iterate(s::MMappedString{ASCII}, i::Int=firstindex(s))
    i > ncodeunits(s) && return nothing
    b = codeunit(s, i)
    Char(b), i + 1
end
Base.@propagate_inbounds function Base.nextind(s::MMappedString{ASCII}, i::Int)
    i == 0 && return 1
    n = ncodeunits(s)
    @boundscheck Base.between(i, 1, n) || throw(Base.BoundsError(s, i))
    i + 1
end
Base.@propagate_inbounds function Base.thisind(s::MMappedString{ASCII}, i::Int)
    i == 0 && return 0
    n = ncodeunits(s)
    i == n + 1 && return i
    @boundscheck Base.between(i, 1, n) || throw(Base.BoundsError(s, i))
    i
end
Base.@propagate_inbounds function Base.prevind(s::MMappedString{ASCII}, i::Int)
    i == 1 && return 0
    n = ncodeunits(s)
    @boundscheck Base.between(i, 1, n + 1) || throw(Base.BoundsError(s, i))
    i - 1
end
@inline function Base.length(s::MMappedString{ASCII}, i::Int, j::Int)
    @boundscheck begin
        0 < i ≤ ncodeunits(s) + 1 || throw(Base.BoundsError(s, i))
        0 ≤ j < ncodeunits(s) + 1 || throw(Base.BoundsError(s, j))
    end
    j < i && return 0
    return j - i + 1
end
Base.isvalid(s::MMappedString{ASCII}, i::Int) = checkbounds(Bool, s, i)
Base.length(s::MMappedString{ASCII}) = length(s.ptr)
Base.lastindex(s::MMappedString{ASCII}) = lastindex(s.ptr)
Base.@propagate_inbounds function Base.getindex(s::MMappedString{ASCII}, i::Int)
    b = codeunit(s, i)
    Char(b)
end
Base.isvalid(s::MMappedString{ASCII}) = Base.isascii(s)
# @Base.propagate_inbounds function Base.getindex(s::MMappedString{ASCII}, r::UnitRange{Int})
#     isempty(r) && return ""
#     return (@view s.ptr[r]) .|> Char |> String
# end
