SuffixArray = Vector{Int32}
DNAString = String

struct GenomeWithSAs
    id::String
    sequence::String
    forwardSA::SuffixArray
    reverseSA::SuffixArray
end

# define equality between GenomeWithSAs
Base.:(==)(x::GenomeWithSAs, y::GenomeWithSAs) = begin
    return x.id == y.id && x.sequence == y.sequence && 
        x.forwardSA == y.forwardSA && x.reverseSA == y.reverseSA;
end

function makeSuffixArray(source::DNAString, circular::Bool)::SuffixArray

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

function makeSuffixArrayT(seqloop::DNAString)::SuffixArray # assumes seqloop is circular

	last::Int32 = trunc(Int32, cld((length(seqloop) + 1) / 2, 3))
	suffixes = Array{SubString}(undef, last * 3)

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

function makeSuffixArray(suffixes::Array{SubString})::SuffixArray

    suffixArray = SuffixArray(undef, length(suffixes))
    suffixArray = sortperm!(suffixArray, suffixes)

    return suffixArray

end

function makeSuffixArrayRanksArray(SA::SuffixArray)::SuffixArray
    len = length(SA)
    RA = SuffixArray(undef, len)
    for i = 1:len
        RA[SA[i]] = i
    end
    return RA
end

import JLD

function writeGenomeWithSAs(filename::String, genome::GenomeWithSAs)
    JLD.jldopen(filename, "w") do file
        write(file, genome.id, genome)
    end
end

function readGenomeWithSAs(filename::String, id::String)::GenomeWithSAs
    JLD.jldopen(filename, "r") do file
        return read(file, id)
    end
end
