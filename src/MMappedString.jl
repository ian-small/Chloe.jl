
module MappedString
export Unicode, ASCII, MMappedString
## Abstract String backed by a readonly Memory Mapped Vector{UInt8}
# as I can't be sure that String(Vector{UInt8}) doesn't make a copy.... :(

# MMappedString simply holds a pointer to a Vector{UInt8} object
# (which is mapped to a mmapped file). All the nasty utf8 logic
# is squirrelled away in the UInt8Utf8 module. We only implement
# the "AbstractString" API here.

using ..UInt8Utf8

import Base
import Base: @propagate_inbounds

struct Unicode end
struct ASCII end

struct MMappedString{T <: Union{Unicode,ASCII}} <: AbstractString
    ptr::Vector{UInt8}
end

MMappedString(ptr) = MMappedString{Unicode}(ptr)

MMappedString(s::String) = MMappedString(Vector{UInt8}(s)) # copy unfortunately

# Must be ASCII! call isvalid(s) to check
MMappedString{ASCII}(s::String) = MMappedString{ASCII}(Vector{UInt8}(s))


Base.:(==)(x::MMappedString, y::MMappedString) = x.ptr == y.ptr

Base.:(==)(x::String, y::MMappedString) = begin
    sa = sizeof(x)
    sa == length(y.ptr) && Base._memcmp(x, pointer(y.ptr), sa) == 0
end

Base.:(==)(y::MMappedString, x::String) = x == y

Base.cmp(a::MMappedString, b::MMappedString) = Base.cmp(a.ptr, b.ptr)

Base.cmp(a::String, b::MMappedString) = begin
    al, bl = sizeof(a), length(b)
    c = Base._memcmp(a, pointer(b.ptr), min(al, bl))
    c < 0 ? -1 : c > 0 ? +1 : Base.cmp(al, bl)
end

@inline Base.cmp(a::MMappedString, b::String) = -Base.cmp(b, a)

@inline Base.hash(s::MMappedString, h::UInt) = Base.hash(s.ptr, h)

Base.show(io::IO, s::MMappedString) = begin
    print(io, "MMappedString[$(length(s.ptr))] @ $(pointer(s.ptr)) of type $(Base.codeunit(s))")
end


Base.pointer(s::MMappedString) = pointer(s.ptr)

Base.pointer(s::MMappedString, i::Integer) = pointer(s.ptr, i)

@propagate_inbounds Base.ncodeunits(s::MMappedString) = utf8_ncodeunits(s.ptr)

@propagate_inbounds Base.codeunit(s::MMappedString) = utf8_codeunit(s.ptr)

@propagate_inbounds Base.codeunit(s::MMappedString, i::Int) = utf8_codeunit(s.ptr, i)

@propagate_inbounds Base.thisind(s::MMappedString, i::Int) = utf8_thisind(s.ptr, i)

@propagate_inbounds Base.nextind(s::MMappedString, i::Int) = utf8_nextind(s.ptr, i)

@propagate_inbounds Base.isvalid(s::MMappedString) = utf8_isvalid(s.ptr)

## required core functionality ##

@propagate_inbounds Base.iterate(s::MMappedString, i::Int=firstindex(s)) = utf8_iterate(s.ptr, i)

@propagate_inbounds Base.getindex(s::MMappedString, i::Int) = utf8_getindex(s.ptr, i)

@propagate_inbounds Base.getindex(s::MMappedString, r::UnitRange{<:Integer}) = s[Int(first(r)):Int(last(r))]

@propagate_inbounds Base.getindex(s::MMappedString, r::UnitRange{Int}) = utf8_getindex(s.ptr, r)

@propagate_inbounds Base.length(s::MMappedString) = utf8_length(s.ptr)

@propagate_inbounds Base.length(s::MMappedString, i::Int, j::Int) = utf8_length(s.ptr, i, j)

## overload methods for efficiency ##

@propagate_inbounds Base.isvalid(s::MMappedString, i::Int) = utf8_isvalid(s.ptr, i)

@propagate_inbounds Base.isascii(s::MMappedString) = utf8_isascii(s.ptr)


# specializations

## Must be ASCII! (Not e.g. latin1). Otherwise the bitwise equality above
# would not work and latin1 from 0x80 up are 2 byte representations in utf-8

Base.show(io::IO, s::MMappedString{ASCII}) = begin
    print(io, "MMappedString[$(length(s.ptr))] @ $(pointer(s.ptr)) of type $(Base.codeunit(s)) ASCII only")
end

@propagate_inbounds Base.iterate(s::MMappedString{ASCII}, i::Int=firstindex(s)) = ascii_iterate(s.ptr, i)

@propagate_inbounds Base.nextind(s::MMappedString{ASCII}, i::Int) = ascii_nextind(s.ptr, i)

@propagate_inbounds Base.thisind(s::MMappedString{ASCII}, i::Int) = ascii_thisind(s.ptr, i)

@propagate_inbounds Base.prevind(s::MMappedString{ASCII}, i::Int) = ascii_prevind(s.ptr, i)

@propagate_inbounds Base.length(s::MMappedString{ASCII}, i::Int, j::Int) = ascii_length(s.ptr, i, j)

@propagate_inbounds Base.isvalid(s::MMappedString{ASCII}, i::Int) = ascii_isvalid(s.ptr, i)

@propagate_inbounds Base.length(s::MMappedString{ASCII}) = length(s.ptr)

@inline Base.firstindex(s::MMappedString{ASCII}) = 1

@propagate_inbounds Base.lastindex(s::MMappedString{ASCII}) = length(s.ptr)

@propagate_inbounds Base.getindex(s::MMappedString{ASCII}, i::Int) = ascii_getindex(s.ptr, i)

@propagate_inbounds Base.isvalid(s::MMappedString{ASCII}) = utf8_isascii(s.ptr)

@propagate_inbounds Base.getindex(s::MMappedString{ASCII}, r::UnitRange{Int}) = ascii_getindex(s.ptr, r)

end
