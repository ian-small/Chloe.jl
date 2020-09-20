
import CodecZlib: GzipDecompressorStream, GzipCompressorStream
import Base

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
            @debug metadata
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
function readFasta(f::IO)::Tuple{String,String}
    seqs = Vector{String}()
    header = strip(readline(f))
    if !startswith(header, ">")
        error("expecting '>' as start of fasta header found: \"$(header)\"")
    end
    nc_id = split(header, " ")[1][2:end]
    for (idx, line) in enumerate(eachline(f))
        line = uppercase(strip(line))
        if length(line) === 0
            continue
        end
        # see https://www.bioinformatics.org/sms/iupac.html
        if match(r"^[ACTGRYNXSWKMBDHV.-]+$", line) === nothing
            error("expecting nucleotide sequence found[$(idx + 1)]: \"$(line)\"")
        end
        push!(seqs, line)
    end
    return string(nc_id), join(seqs, "")
end

const COMP = Dict('A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', 
                  'R' => 'Y', 'Y' => 'R', 'N' => 'N', 'X' => 'X')
    
function revComp(dna::AbstractString)::String # where {T <: AbstractString}
    reverse(map(x -> get(COMP, x, 'N'), dna))
end

@inline function frameCounter(base::T, addition::T) where {T <: Integer}
    result = (base - addition) % 3
    if result <= 0
        result = 3 + result
    end
    result
end

@inline function phaseCounter(base::Int8, addition::Int32)::Int8
    result = (base - addition) % 3
    if result < 0
        result = 3 + result
    end
    result
end

@inline function rangesOverlap(start1::T, length1::T, start2::T, length2::T)::Bool where {T <: Integer}
    if start1 >= start2 + length2 || start2 >= start1 + length1
        return false
    end
    true
end


@inline function genome_wrap(genome_length::T, position::T)::T where {T <: Integer}
    if 0 < position <= genome_length
        return position
    end
    while position <= 0
        position = position + genome_length
    end
    while position > genome_length
        position = position - genome_length
    end
    position
end

# wraps range within genome loop
Base.@propagate_inbounds function range_wrap!(range::UnitRange{T}, genome_length::T)::UnitRange{T} where {T <: Integer}
    @boundscheck range.start <= range.stop  && genome_length > 0 || throw(BoundError("start after stop"))
    loop_length = genome_length + genome_length - 1
    # if start of range is negative, move range to end of genome
    while range.start <= 0
        lengthminus1 = range.stop - range.start
        range.start += genome_length
        range.stop = range.start + lengthminus1
    end
    # if end of range is beyond end of loop, move range to beginning of loop
    while range.stop > loop_length
        range.start -= genome_length
        range.stop -= genome_length
    end
    # if start of range is beyond end of genome, move range to beginning of genome
    while range.start > genome_length
        range.start -= genome_length
        range.stop -= genome_length
    end
    @boundscheck 0 < range.start <= genome_length &&  0 < range.stop <= loop_length || throw(BoundsError("out of range"))
    return range
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
