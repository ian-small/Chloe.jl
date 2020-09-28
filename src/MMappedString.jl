
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
# Base.eltype(::Type{MMappedString}) = Char

MMappedString(ptr) = MMappedString{Unicode}(ptr)

MMappedString(s::String) = MMappedString(Vector{UInt8}(s)) # copy unfortunately

# Must be ASCII! call isvalid(s) to check
MMappedString{ASCII}(s::String) = MMappedString{ASCII}(Vector{UInt8}(s))

@inline Base.:(==)(x::MMappedString, y::MMappedString) = x.ptr == y.ptr
@inline Base.cmp(a::MMappedString, b::MMappedString) = Base.cmp(a.ptr, b.ptr)

Base.:(==)(x::String, y::MMappedString) = begin
    sa = sizeof(x)
    sa == length(y.ptr) && Base._memcmp(x, pointer(y.ptr), sa) == 0
end

Base.:(==)(y::MMappedString, x::String) = x == y


Base.cmp(a::String, b::MMappedString) = begin
    al, bl = sizeof(a), length(b.ptr)
    c = Base._memcmp(a, pointer(b.ptr), min(al, bl))
    c < 0 ? -1 : c > 0 ? +1 : Base.cmp(al, bl)
end

@inline Base.cmp(a::MMappedString, b::String) = -Base.cmp(b, a)

@inline Base.hash(s::MMappedString, h::UInt) = Base.hash(s.ptr, h)

Base.show(io::IO, s::MMappedString) = begin
    print(io, "MMappedString[$(length(s.ptr))] @ $(pointer(s.ptr))")
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
# @propagate_inbounds getindex(s::MMappedString, r::UnitRange{Int}) = SubString(s, r)

@propagate_inbounds Base.length(s::MMappedString) = utf8_length(s.ptr)

@propagate_inbounds Base.length(s::MMappedString, i::Int, j::Int) = utf8_length(s.ptr, i, j)

## overload methods for efficiency ##

@propagate_inbounds Base.isvalid(s::MMappedString, i::Int) = utf8_isvalid(s.ptr, i)

@propagate_inbounds Base.isascii(s::MMappedString) = utf8_isascii(s.ptr)

@propagate_inbounds Base.String(s::MMappedString) = begin
    parent = s.ptr
    return GC.@preserve parent Base.unsafe_string(pointer(parent), length(s.ptr))
end

# substrings 

@propagate_inbounds Base.String(s::SubString{MMappedString}) = begin
    parent = s.string.ptr
    return GC.@preserve parent Base.unsafe_string(pointer(parent, s.offset + 1), s.ncodeunits)
end

function Base.cmp(sa::SubString{Union{MMappedString,String}}, sb::Union{String,SubString{MMappedString}})

    na = sizeof(sa)
    nb = sizeof(sb)
    c = Base._memcmp(pointer(sa), pointer(sb), min(na, nb))
    return c < 0 ? -1 : c > 0 ? +1 : cmp(na, nb)
end

# reverse comparisons
Base.cmp(sb::SubString{MMappedString}, sa::SubString{String}) = -Base.cmp(sa, sb)
Base.cmp(sb::String, sa::SubString{MMappedString}) = -Base.cmp(sa, sb)

Base.pointer(x::SubString{MMappedString}) = pointer(x.string.ptr) + x.offset
Base.pointer(x::SubString{MMappedString}, i::Integer) = pointer(x.string.ptr) + x.offset + (i - 1)

# specializations

## Must be ASCII! (Not e.g. latin1). Otherwise the bitwise equality above
# would not work and latin1 from 0x80 up are 2 byte representations in utf-8

# The use of MMappedString{ASCII} is that now length, lastindex, m[i] *won't* need to
# do a scan of the bytes to find the number of characters

Base.show(io::IO, s::MMappedString{ASCII}) = begin
    print(io, "MMappedString[$(length(s.ptr))] @ $(pointer(s.ptr)) ASCII only")
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

# the problem here is that MMappedString{Unicode} might
# actually be ASCII so x.ptr == y.ptr will be *way* quicker than scanning the characters.
# you can have all charaters the same:
# all([a == b for (a,b) in zip(string, mmap)) && length(string) == length(mmap)
# *and* string != mmap
function allow_latin1()
    @eval begin
        @propagate_inbounds Base.:(==)(m::MMappedString{ASCII}, s::String) = latin1_eq(m.ptr, s)
        @propagate_inbounds Base.:(==)(s::String, m::MMappedString{ASCII}) =  m == s
        @propagate_inbounds Base.cmp(m::MMappedString{ASCII}, s::String) = latin1_cmp(m.ptr, s)
        @propagate_inbounds Base.cmp(s::String, m::MMappedString{ASCII}) = -Base.cmp(m, s)
    end
end

## Regex

# subvert the type checks of julia by disguising MMappedString as a SubString
Base.occursin(re::Regex , s::MMappedString; offset::Integer=0) = begin
    Base.occursin(re, SubString(s, 1); offset=offset)
end
function startswith(s::MMappedString, r::Regex)
    Base.startswith(SubString(s, 1), r)
end
Base.findnext(r::Regex, s::MMappedString, idx::Integer) = begin
    Base.findnext(r, SubString(s, 1), idx)
end
# don't know why I have to enumerate subclasses
function Base.match(re::Regex, str::Union{SubString{MMappedString{Unicode}},SubString{MMappedString{ASCII}},MMappedString}, idx::Integer, add_opts::UInt32=UInt32(0))
    Base.compile(re)
    opts = re.match_options | add_opts
    matched, data = Base.PCRE.exec_r_data(re.regex, str, idx - 1, opts)
    if !matched
        Base.PCRE.free_match_data(data)
        return nothing
    end
    n = div(Base.PCRE.ovec_length(data), 2) - 1
    p = Base.PCRE.ovec_ptr(data)
    mat = SubString(str, Base.unsafe_load(p, 1) + 1, Base.prevind(str, Base.unsafe_load(p, 2) + 1))
    # needs to be SubString{String} since this is the type of captures field :(
    cap = Union{Nothing,SubString{String}}[Base.unsafe_load(p, 2i + 1) == Base.PCRE.UNSET ? nothing :
                                        SubString(str, Base.unsafe_load(p, 2i + 1) + 1,
                                                Base.prevind(str, Base.unsafe_load(p, 2i + 2) + 1)) for i = 1:n]
    off = Int[ Base.unsafe_load(p, 2i + 1) + 1 for i = 1:n ]
    result = RegexMatch(mat, cap, Base.unsafe_load(p, 1) + 1, off, re)
    Base.PCRE.free_match_data(data)
    return result
end

function Base.PCRE.exec(re, subject::Union{SubString{MMappedString{Unicode}},SubString{MMappedString{ASCII}},MMappedString}, offset, options, match_data)
    rc = ccall((:pcre2_match_8, Base.PCRE.PCRE_LIB), Cint,
               (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Csize_t, UInt32, Ptr{Cvoid}, Ptr{Cvoid}),
               re, subject, ncodeunits(subject), offset, options, match_data, Base.PCRE.get_local_match_context())
    # rc == -1 means no match, -2 means partial match.
    rc < -2 && error("PCRE.exec error: $(err_message(rc))")
    return rc >= 0
end

end # module