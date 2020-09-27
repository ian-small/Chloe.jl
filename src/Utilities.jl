
import CodecZlib: GzipDecompressorStream, GzipCompressorStream
import Base
struct FwdRev{T}
    forward::T
    reverse::T
end

Base.:(==)(x::FwdRev{T}, y::FwdRev{T}) where T = x.forward == y.forward && x.reverse == y.reverse


datasize(t::T) where T = sizeof(t)
datasize(f::FwdRev{T}) where T = sizeof(FwdRev{T}) + datasize(f.forward) + datasize(f.reverse)
datasize(v::Vector{T}) where T = length(v) == 0 ? 0 : sum(datasize(a) for a in v)
datasize(t::Dict{K,V}) where {K,V} = length(t) == 0 ? 0 : sum(datasize(e.first) + datasize(e.second) for e in t)

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

function readFasta(fasta::String)::Tuple{String,String}
    if !isfile(fasta)
        error("$(fasta): not a file!")
    end
    maybe_gzread(fasta) do io
        readFasta(io)
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


@inline function rangesOverlap(start1::Int32, length1::Int32, start2::Int32, length2::Int32)::Bool
    if start1 >= start2 + length2 || start2 >= start1 + length1
        return false
    end
    true
end

function isStartCodon(codon::AbstractString, allow_editing::Bool, allow_GTG::Bool)::Bool
    codon == "ATG" && return true
    if allow_editing && codon == "ACG"
        return true
    elseif allow_GTG && codon == "GTG"
        return true
    end
    return false
end

function isStopCodon(codon::AbstractString, allow_editing::Bool)::Bool
    if codon in ["TAA","TAG","TGA"]
        return true
    elseif allow_editing && codon in ["CAA","CAG","CGA"]
        return true
    end
    return false
end

const GENETIC_CODE = Dict{String,Char}(
    "TTT" => 'F',"TTC" => 'F',"TTA" => 'L',"TTG" => 'L',
    "CTT" => 'L',"CTC" => 'L',"CTA" => 'L',"CTG" => 'L',
    "ATT" => 'I',"ATC" => 'I',"ATA" => 'I',"ATG" => 'M',
    "GTT" => 'V',"GTC" => 'V',"GTA" => 'V',"GTG" => 'V',
    "TCT" => 'S',"TCC" => 'S',"TCA" => 'S',"TCG" => 'S',
    "CCT" => 'P',"CCC" => 'P',"CCA" => 'P',"CCG" => 'P',
    "ACT" => 'T',"ACC" => 'T',"ACA" => 'T',"ACG" => 'T',
    "GCT" => 'A',"GCC" => 'A',"GCA" => 'A',"GCG" => 'A',
    "TAT" => 'Y',"TAC" => 'Y',"TAA" => '*',"TAG" => '*',
    "CAT" => 'H',"CAC" => 'H',"CAA" => 'Q',"CAG" => 'Q',
    "AAT" => 'N',"AAC" => 'N',"AAA" => 'K',"AAG" => 'K',
    "GAT" => 'D',"GAC" => 'D',"GAA" => 'E',"GAG" => 'E',
    "TGT" => 'C',"TGC" => 'C',"TGA" => '*',"TGG" => 'W',
    "CGT" => 'R',"CGC" => 'R',"CGA" => 'R',"CGG" => 'R',
    "AGT" => 'S',"AGC" => 'S',"AGA" => 'R',"AGG" => 'R',
    "GGT" => 'G',"GGC" => 'G',"GGA" => 'G',"GGG" => 'G')

function translateDNA(dna::AbstractString)::String
    @inbounds String([get(GENETIC_CODE, SubString(dna, i, i + 2), 'X') for i in 1:3:length(dna) - 2])
end
function translateDNA(dna::AbstractString, start::Integer, stop::Integer)::String
    String([get(GENETIC_CODE, SubString(dna, i, i + 2), 'X') for i in start:3:stop])
end

const COMP = Dict('A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', 
                  'R' => 'Y', 'Y' => 'R', 'N' => 'N', 'X' => 'X')
    
function revComp(dna::AbstractString)::String
    String(reverse(map(x -> get(COMP, x, 'N'), dna)))
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

struct IFasta
    io::IO
    name::String
    full::Bool # full header
    close::Bool
end

Base.eltype(::Type{IFasta}) = Tuple{String,String}
Base.IteratorSize(::Type{IFasta}) = Base.SizeUnknown()

function iterFasta(fasta::String; full::Bool=false)
    if !isfile(fasta)
        error("$(fasta): not a file!")
    end
    io = if endswith(fasta, ".gz")
        open(fasta) |> GzipDecompressorStream
    else
        open(fasta)
    end
    return IFasta(io, fasta, full, true)
end 

iterFasta(io::IO, name::String="<stream>", close::Bool=true; full::Bool=false) = IFasta(io, name, full, close)


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

function Base.iterate(i::IFasta, state=nothing)
    done = () -> i.close && close(i.io)

    if state === nothing
        if eof(i.io)
            done()
            return nothing
        end
        state = findStart(i.io)
        if state == "" # a single line of spaces is i guess valid
            done()
            return nothing
        end
        lineno, header = 1, state
    else
        lineno, header = state
        if header === nothing
            done()
            return nothing
        end
    end

    if !startswith(header, ">")
        done()
        error("$(i.name): expecting \">\" at line $(lineno) got: $(str_truncate(header))")
    end

    id = String(i.full ? header[2:end] : split(header, " ")[1][2:end])
    
    if length(id) == 0
        done()
        error("$(i.name): no FASTA ID at line $(lineno)")
    end

    seq, lineno, more = readFastaBody(i, lineno)
    return ((id, seq), (lineno, more))
end

function readFastaBody(f::IFasta, lineno::Int=0)::Tuple{String,Int,Union{String,Nothing}}
    seqs = Vector{String}()
    while !eof(f.io)
        line = strip(readline(f.io))
        lineno += 1
        if length(line) === 0
            continue
        end
        if startswith(line, ">")
            return join(seqs, ""), lineno, line
        end
        line = uppercase(line)
        # see https://www.bioinformatics.org/sms/iupac.html
        if match(r"^[ACTGRYNXSWKMBDHV.-]+$", line) === nothing
            error("$(f.name): expecting nucleotide sequence found[$(lineno)]: \"$(str_truncate(line))\"")
        end
        push!(seqs, line)
    end
    return join(seqs, ""), lineno, nothing
end
