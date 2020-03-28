SuffixVector = Vector{Int32}
SuffixArray = Array{Int32}

struct GenomeWithSAs
    id::String
    sequence::String
    forwardSA::SuffixVector
    reverseSA::SuffixVector
end

function makeSuffixArray(source, circular)::SuffixArray

    if circular
		last = Int((length(source) + 1) / 2)
    else
        last = length(source)
    end

    suffixes = Array{SubString}(undef, last)
    for offset = 1:last
        suffixes[offset] = SubString(source, offset)
    end

	suffixArray = SuffixArray(undef, last)
    suffixArray = sortperm!(suffixArray, suffixes)

    return suffixArray

end

function makeSuffixArrayT(seqloop::AbstractString)::SuffixArray # assumes seqloop is circular

	last::Int = trunc(Int32, cld((length(seqloop) + 1) / 2, 3))
	suffixes = Vector{SubString}(undef, last * 3)

	frame = translateDNA(seqloop)
	for offset = 1:last
		suffixes[offset] = SubString(frame, offset)
	end
	seqloop = seqloop[2:end] * seqloop[1]
	frame = translateDNA(seqloop)
	for offset = last + 1:last * 2
		suffixes[offset] = SubString(frame, offset)
	end
	seqloop = seqloop[2:end] * seqloop[1]
	frame = translateDNA(seqloop)
	for offset = last * 2 + 1:last * 3
		suffixes[offset] = SubString(frame, offset)
	end
	return makeSuffixArray(suffixes)
end

function makeSuffixArray(suffixes::Vector{SubString})::SuffixArray

    suffixArray = SuffixArray(undef, length(suffixes))
    suffixArray = sortperm!(suffixArray, suffixes)

    return suffixArray

end

function makeSuffixArrayRanksArray(SA::Union{SuffixVector, SuffixVector})::Array{Int32}
    len = length(SA)
    RA = Array{Int32}(undef, len)
    for i = 1:len
        RA[SA[i]] = i
    end
    return RA
end

using JLD

function writeGenomeWithSAs(filename::String, genome::GenomeWithSAs)
    jldopen(filename, "w") do file
        write(file, genome.id, genome)
    end
end

function readGenomeWithSAs(filename::String, id::String)::GenomeWithSAs
    jldopen(filename, "r") do file
        return read(file, id)
    end
end
