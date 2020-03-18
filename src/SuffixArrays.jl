struct GenomeWithSAs
    id::String
    sequence::String
    forwardSA::Vector{Int32}
    reverseSA::Vector{Int32}
end

function makeSuffixArray(source, circular)

    if circular
		last = Int((length(source) + 1) / 2)
    else
        last = length(source)
    end

    suffixes = Array{SubString}(undef, last)
    for offset = 1:last
        suffixes[offset] = SubString(source, offset)
    end

	suffixArray = Array{Int32}(undef, last)
    suffixArray = sortperm!(suffixArray, suffixes)

    return suffixArray

end

function makeSuffixArrayT(seqloop) # assumes seqloop is circular

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

function makeSuffixArray(suffixes::Vector{SubString})

	suffixArray = Array{Int32}(undef, length(suffixes))
    suffixArray = sortperm!(suffixArray, suffixes)

    return suffixArray

end

function makeSuffixArrayRanksArray(SA)
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

function readGenomeWithSAs(filename::String, id::String)
    jldopen(filename, "r") do file
        return read(file, id)
    end
end
