
module MappedString
export Unicode, ASCII, MMappedString
## Abstract String backed by a readonly Memory Mapped Vector{UInt8}
# as I can't be sure that String(Vector{UInt8}) doesn't make a copy.... :(

using ..UInt8Utf8


struct Unicode end
struct ASCII end
struct MMappedString{T <: Union{Unicode,ASCII}} <: AbstractString
    ptr::Vector{UInt8}
end

MMappedString(ptr) = MMappedString{Unicode}(ptr)

MMappedString(s::String) = MMappedString(Vector{UInt8}(s)) # copy unfortunately
MMappedString{ASCII}(s::String) = MMappedString{ASCII}(Vector{UInt8}(s))

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
    print(io, "MMappedString[$(length(s.ptr))] @ $(pointer(s.ptr)) of type $(Base.codeunit(s))")
end

# Base.print(io::IO, s::MMappedString) = print(io, string(s))
# Base.textwidth(s::MMappedString) = textwidth(string(s))


# Base.string(x::MMappedString) = String(x)

# Base.convert(::Type{MMappedString}, x::String) = MMappedString(x)
# Base.convert(::Type{String}, x::MMappedString) = String(x)
# Base.String(x::MMappedString) = String(copy(x.ptr))

Base.pointer(s::MMappedString) = pointer(s.ptr)
Base.pointer(s::MMappedString, i::Integer) = pointer(s.ptr, i)

Base.ncodeunits(s::MMappedString) = utf8_ncodeunits(s.ptr)
Base.codeunit(s::MMappedString) = utf8_codeunit(s.ptr)
Base.codeunit(s::MMappedString, i::Int) = utf8_codeunit(s.ptr, i)


Base.thisind(s::MMappedString, i::Int) = utf8_thisind(s.ptr, i)

Base.nextind(s::MMappedString, i::Int) = utf8_nextind(s.ptr, i)

Base.isvalid(s::MMappedString) = utf8_isvalid(s.ptr)

## required core functionality ##

Base.iterate(s::MMappedString, i::Int=firstindex(s)) = utf8_iterate(s.ptr, i)


Base.getindex(s::MMappedString, i::Int) = utf8_getindex(s.ptr, i)


Base.getindex(s::MMappedString, r::UnitRange{<:Integer}) = s[Int(first(r)):Int(last(r))]

Base.getindex(s::MMappedString, r::UnitRange{Int}) = utf8_getindex(s.ptr, r)


Base.length(s::MMappedString) = utf8_length(s.ptr)

Base.length(s::MMappedString, i::Int, j::Int) = utf8_length(s.ptr, i, j)


## overload methods for efficiency ##

Base.isvalid(s::MMappedString, i::Int) = utf8_isvalid(s.ptr, i)

Base.isascii(s::MMappedString) = utf8_isascii(s.ptr)


# specializations

## Must be ASCII! Otherwise the bitwise equality above
# would not work and latin1 from 0x80 up are 2 byte representations in utf-8


function Base.show(io::IO, s::MMappedString{ASCII})
    print(io, "MMappedString[$(length(s.ptr))] @ $(pointer(s.ptr)) of type $(Base.codeunit(s)) Latin1")
end
Base.iterate(s::MMappedString{ASCII}, i::Int=firstindex(s)) = ascii_iterate(s.ptr, i)

# step forward and back by one byte == one character
Base.nextind(s::MMappedString{ASCII}, i::Int) = ascii_nextind(s.ptr, i)

Base.thisind(s::MMappedString{ASCII}, i::Int) = ascii_thisind(s.ptr, i)

Base.prevind(s::MMappedString{ASCII}, i::Int) = ascii_prevind(s.ptr, i)


Base.length(s::MMappedString{ASCII}, i::Int, j::Int) = ascii_length(s.ptr, i, j)

Base.isvalid(s::MMappedString{ASCII}, i::Int) = ascii_isvalid(s.ptr, i)
Base.length(s::MMappedString{ASCII}) = length(s.ptr)
Base.firstindex(s::MMappedString{ASCII}) = 1
Base.lastindex(s::MMappedString{ASCII}) = length(s.ptr)
Base.getindex(s::MMappedString{ASCII}, i::Int) = ascii_getindex(s.ptr, i)

Base.isvalid(s::MMappedString{ASCII}) = utf8_isascii(s.ptr)

end
