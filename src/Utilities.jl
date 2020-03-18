
function gbff2fasta(infile)
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
            println(metadata)
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
    return
end

function readFasta(file)
    id = ""
    seq = ""
    open(file) do f
        header = readline(f)
        id = split(header, " ")[1][2:end]
        while !eof(f)
            seq = seq * uppercase(readline(f))
        end
    end
    return id, seq
end

function revComp(dna)
    comp = Dict('A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', 'R' => 'Y', 'Y' => 'R', 'N' => 'N', 'X' => 'X')
    reverse(map(x->get(comp, x, 'N'), dna))
end

function frameCounter(base::Integer, addition::Integer)
    result = (base - addition) % 3
    if result <= 0
        result = 3 + result
    end
    return result
end

function phaseCounter(base::Integer, addition::Integer)
    result = (base - addition) % 3
    if result < 0
        result = 3 + result
    end
    return result
end

function rangesOverlap(start1::Integer, length1::Integer, start2::Integer, length2::Integer)
    if start1 >= start2 + length2 || start2 >= start1 + length1
        return false
    else
        return true
    end
end

# wraps to genome length
function genome_wrap(genome_length::Integer, position::Integer)
    0 < position <= genome_length && return position
    position <= 0 && return genome_length + position
    position > genome_length && return position - genome_length
end

# wraps to loop length, i.e. allows position to exceed genome length
function loop_wrap(genome_length::Integer, position::Integer)
    0 < position <= genome_length + genome_length - 1 && return position
    position <= 0 && return genome_length + position
    return position - genome_length
end

# wraps range within genome loop
function range_wrap(genome_length::Integer, range::UnitRange{Integer})
    @assert range.start <= range.stop
    loop_length = genome_length + genome_length - 1
    # if start of range is negative, move range to end of genome
    if range.start <= 0
        lengthminus1 = range.stop - range.start
        range.start = genome_length + range.start
        range.stop = range.start + lengthminus1
    end
    # if end of range is beyond end of loop, move range to beginning of loop
    if range.stop > loop_length
        range.start = range.start - genome_length
        range.stop = range.stop - genome_length
    end
    # if start of range is beyond end of genome, move range to beginning of genome
    if range.start > genome_length
        range.start = range.start - genome_length
        range.stop = range.stop - genome_length
    end
    @assert 0 < range.start <= genome_length
    @assert 0 < range.stop <= loop_length
    return range::UnitRange{Integer}
end

function isStartCodon(codon::AbstractString, allow_editing::Bool, allow_GTG::Bool)
    codon == "ATG" && return true
    if allow_editing && codon == "ACG"
        return true
    elseif allow_GTG && codon == "GTG"
        return true
    end
    return false
end

function isStopCodon(codon::AbstractString, allow_editing::Bool)
    if codon in ["TAA","TAG","TGA"]
        return true
    elseif allow_editing && codon in ["CAA","CAG","CGA"]
        return true
    end
    return false
end

const genetic_code = Dict{String,Char}("TTT" => 'F',"TTC" => 'F',
    "TTA" => 'L',"TTG" => 'L',
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

function translateDNA(dna::AbstractString)
    peptide_length = fld(length(dna), 3)
    peptide = Array{Char}(undef, peptide_length)
    aa = 0
    for i = 1:3:length(dna) - 2
        aa += 1
        peptide[aa] = get(genetic_code, SubString(dna, i, i + 2), 'X')
    end
    return String(peptide)
end
