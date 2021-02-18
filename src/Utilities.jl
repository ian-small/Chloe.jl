import CodecZlib: GzipDecompressorStream, GzipCompressorStream
import Base

import Printf: @sprintf

const REENTRANT_LOCK = ReentrantLock()

struct FwdRev{T}
    forward::T
    reverse::T
end

Base.:(==)(x::FwdRev{T}, y::FwdRev{T}) where T = x.forward == y.forward && x.reverse == y.reverse

datasize(t::T) where T = sizeof(t)
datasize(f::FwdRev{T}) where T = sizeof(FwdRev{T}) + datasize(f.forward) + datasize(f.reverse)
datasize(v::Vector{T}) where T = length(v) == 0 ? 0 : sum(datasize(a) for a in v)
datasize(t::Dict{K,V}) where {K,V} = length(t) == 0 ? 0 : sum(datasize(e.first) + datasize(e.second) for e in t)

struct CircularVector
    v::Vector{Int32}
end

@inline Base.length(cv::CircularVector) = Int32(length(cv.v))
@inline Base.getindex(cv::CircularVector, i::Int32) = @inbounds getindex(cv.v, mod1(i, length(cv)))
function Base.getindex(cv::CircularVector, r::UnitRange{<:Integer})
    len = length(cv)
	start = mod1(r.start, len)
    stop = mod1(r.stop, len)
	if start < stop
		return cv.v[start:stop]
    else
        return vcat(cv.v[start:end],cv.v[1:stop])
    end
end
@inline Base.setindex!(cv::CircularVector, value::Int32, i::Int32) = @inbounds setindex!(cv.v, value, mod1(i, length(cv)))
@inline Base.push!(cv::CircularVector, value::Int32) = @inbounds push!(cv.v, value)
function Base.sum(cv::CircularVector, r::UnitRange{<:Integer})
    sum = 0
    for i in r
        sum += cv[i]
    end
    return sum
end

using BioSequences

struct CircularSequence
    sequence::LongDNASeq
end

@inline Base.length(cs::CircularSequence) = Int32(length(cs.sequence))
@inline Base.getindex(cs::CircularSequence, i::Int32) = @inbounds getindex(cs.sequence, mod1(i, length(cs)))

function Base.getindex(cs::CircularSequence, r::UnitRange{<:Integer})
    len = length(cs.sequence)
	start = mod1(r.start, len)
    stop = mod1(r.stop, len)
	if start < stop
		return cs.sequence[start:stop]
    else
        return append!(cs.sequence[start:end],cs.sequence[1:stop])
    end
end

function reverse_complement(cs::CircularSequence)
    return CircularSequence(BioSequences.reverse_complement(cs.sequence))
end

@inline function getcodon(cs::CircularSequence, index::Int32)
    return (cs[index], cs[index + Int32(1)], cs[index + Int32(2)])
end

function gbff2fasta(infile::String)
    open(infile) do f
        while !eof(f)
            line = readline(f)
            fields = split(line)
            accession = fields[2]
            metadata = fields[2] * "\t" * fields[3] * "\t"
            while !occursin("ORGANISM", line)
                line = readline(f)
            end
            metadata = metadata * line[13:end] * "\t"
            line = readline(f)
            while !occursin("REFERENCE", line)
                metadata = metadata * line[13:end]
                line = readline(f)
            end
            while !startswith(line, "ORIGIN")
                line = readline(f)
            end
            open(accession * ".fna", "w") do o
                write(o, ">", accession, "\n")
                while !startswith(line, "//")
                    line = readline(f)
                    write(o, uppercase(join(split(line[11:end]))))
                end
                write(o, "\n")
            end
        end
    end
end

# const ns(td) = Time(Nanosecond(td))
ns(td) = @sprintf("%.3fs", td / 1e9)
elapsed(st) = @sprintf("%.3fs", (time_ns() - st) / 1e9)

function human(num::Integer)::String
    if num == 0
        return "0B"
    end
    magnitude = floor(Int, log10(abs(num)) / 3)
    val = num / (1000^magnitude)
    sval = @sprintf("%.1f", val)
    if magnitude > 7
        return "$(sval)YB"
    end
    p = ["", "k", "M", "G", "T", "P", "E", "Z"][magnitude + 1]
    return "$(sval)$(p)B"
end

function maybe_gzread(f::Function, filename::String)
    if endswith(filename, ".gz")
        open(z -> z |> GzipDecompressorStream |> f, filename)
    else
        open(f, filename)
    end
end

function maybe_gzwrite(f::Function, filename::String)
    
    function gzcompress(f::Function, fp::IO)
        o = GzipCompressorStream(fp)
        try
            f(o)
        finally
            close(o)
        end
    end

    if endswith(filename, ".gz")
        open(fp -> gzcompress(f, fp), filename, "w")
    else
        open(f, filename, "w")
    end
end


function readFasta(f::IO, name::String="<stream>")::Tuple{String,String}
    for res in iterFasta(f, name)
        return res
    end

    error("$(name): no FASTA data!")
end

function readFasta(fasta::String)::Tuple{String,String}
    for res in iterFasta(fasta)
        return res
    end

    error("$(fasta): no FASTA data!")

end


@inline function frameCounter(base::Int8, addition::Int32)
    result = (base - addition) % 3
    if result <= 0
        result = 3 + result
    end
    result
end

@inline function phaseCounter(base::Int8, addition::Integer)::Int8
    result::Int8 = (base - addition) % 3
    if result < 0
        result = Int8(3) + result
    end
    result
end

@inline function genome_wrap(genome_length::Integer, position::Integer)
    if 0 < position ≤ genome_length
    return position
    end
    while position > genome_length
        position -= genome_length
    end 
    while position ≤ 0
        position += genome_length
    end
    position
end

@inline function circulardistance(start, stop, seqlength)
    return stop ≥ start ? stop - start : stop + seqlength - start
end


@inline function rangesOverlap(start1::Int32, length1::Int32, start2::Int32, length2::Int32)::Bool
    if start1 >= start2 + length2 || start2 >= start1 + length1
        return false
    end
    true
end


# --- wrapped iterators

struct ModuloIteratorState{T <: Integer}
    start::T
    step::T
    length::T # number of iterations
    gene_length::T
end
struct VModuloIteratorState{T <: Integer,V}
    start::T
    step::T
    length::T # number of iterations
    v::Vector{V}
end

Base.length(t::ModuloIteratorState{T}) where {T} = t.length
Base.eltype(::Type{ModuloIteratorState{T}}) where {T} = T
Base.IteratorSize(::Type{ModuloIteratorState{T}}) where {T} = Base.HasLength()

Base.length(t::VModuloIteratorState{T,V}) where {T,V} = t.length
Base.eltype(::Type{VModuloIteratorState{T,V}}) where {T,V} = V
Base.IteratorSize(::Type{VModuloIteratorState{T,V}}) where {T,V} = Base.HasLength()

function Base.iterate(i::ModuloIteratorState{T}, state=(i.start, i.length)) where {T}
    pos, n = state
        n -= 1
    if n < 0
        return nothing
    end
    pos += i.step
    m = i.gene_length
    if pos > m
        pos -= m
    elseif pos ≤ 0
        pos += m
    end
    return (pos, (pos, n))
end
function Base.iterate(i::VModuloIteratorState{T,V}, state=(i.start, i.length)) where {T,V}
    pos, n = state
        n -= 1
    if n < 0
        return nothing
    end
    pos += i.step
    m = length(i.v)
    if pos > m
        pos -= m
    elseif pos ≤ 0
        pos += m
    end
    @inbounds v = i.v[pos]
    return (v, (pos, n))
end
function iter_wrap(r::UnitRange{<:Integer}, gene_length::T) where {T <: Integer}
    # for i in iter_wrap(-5:5, 10)
    # 5,6,7,8,9,10,1,2,3,4,5
    # end
    @boundscheck 1 ≤ gene_length || throw(Base.BoundsError("step size too large"))
    start = genome_wrap(gene_length, r.start)
    return ModuloIteratorState(promote(start - 1, 1, length(r), gene_length)...)
end
function iter_wrap(r::UnitRange{<:Integer}, v::Vector{T}) where {T}
    @boundscheck 1 ≤ length(v) || throw(Base.BoundsError("step size too large"))
    start = genome_wrap(length(v), r.start)
    return VModuloIteratorState(start - 1, 1, length(r), v)
end
function iter_wrap(r::StepRange{<:Integer,<:Integer}, gene_length::T) where {T <: Integer}
    @boundscheck abs(r.step) ≤ gene_length || throw(Base.BoundsError("step size too large"))
    return ModuloIteratorState(promote(genome_wrap(gene_length, r.start) - r.step, r.step, length(r), gene_length)...)
end
function iter_wrap(r::StepRange{<:Integer,<:Integer}, v::Vector{T}) where {T}
    @boundscheck abs(r.step) ≤ length(v) || throw(Base.BoundsError("step size too large"))
    start = genome_wrap(length(v), r.start)
    return VModuloIteratorState(start - r.step, r.step, length(r), v)
end

function findStart(io::IO)
    while !eof(io)
        h = strip(readline(io))
        if length(h) > 0
            return h
        end
    end
    ""
end

function str_truncate(s::String, width=80)
    n = length(s)
    if n <= width
        s
    else
        s[1:thisind(s, width)] * "..."
    end
end

function findfastafile(dir::String, base::AbstractString)::Union{String,Nothing}
    suffixes = [".fa", ".fna", ".fasta"]
    for suffix in suffixes
        path = normpath(joinpath(dir, base * suffix))
        if isfile(path)
            return path
        end
    end
    return nothing
end
