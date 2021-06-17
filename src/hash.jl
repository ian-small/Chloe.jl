using BioSequences
using FASTX
using Statistics

const KMERSIZE = 10
const SKETCHSIZE = 400

function minhash_references(; fasta_files=Vector{String}, output=output)
    seqs = Vector{FASTA.Record}()
    for file in fasta_files
        println(file)
        if isdir(file)
            minhash_references(fasta_files=filter(x -> endswith(x, ".fasta"), readdir(file, join=true)), output=output)
        end
        reader = open(FASTA.Reader, file)
        record = FASTA.Record() 
        read!(reader, record)
        push!(seqs, record)
        close(reader)
    end

    hashes = Dict{String,Vector{Int64}}()
    for seqrec in seqs
        seq = FASTA.sequence(seqrec)
        # calculate entropy mask
        emask = entropy_mask(CircularSequence(seq), Int32(KMERSIZE))
        # convert low entropy regions to Ns
        for (i, bit) in enumerate(emask)
            if bit; seq[i] = DNA_N; end
        end
        # hash
        hashes[FASTA.identifier(seqrec)] = minhash(seq, KMERSIZE, SKETCHSIZE).sketch
    end
    writeminhashes(output, hashes)
end

function writeminhashes(filename::String, hashes::Dict{String,Vector{Int64}})
    outfile = open(filename, "w")
    for (id, hash) in hashes
        write(outfile, UInt8(length(id)), id, hash)
    end
    close(outfile)
end

function readminhashes(infile::String)::Dict{String,Vector{Int64}}
    hashes = Dict{String,Vector{Int64}}()
    # add sizehint! once number of references is established
    open(infile, "r") do hashfile
        while !eof(hashfile)
            length_id = read(hashfile, UInt8)
            id = String(read(hashfile, length_id, all=true))
            h = Vector{Int64}(undef, SKETCHSIZE)
            hashes[id] = read!(hashfile, h)
        end
    end
    return hashes
end

function searchhashes(hash, hashes)::Vector{Tuple{String,Int}}
    scores = Vector{Tuple{String,Int}}(undef, length(hashes))
    for (index, (id, h)) in enumerate(hashes)
        scores[index] = (id, hashcount(hash, h))
    end
    return sort!(scores, by=last, rev=true)
end

function hashcount(hash1, hash2)::Int
    l = length(hash1)
    i = 1
    j = 1
    count = 0
    while l > i && l > j
        if hash1[i] < hash2[j]
            i += 1
        elseif hash1[i] > hash2[j]
            j += 1
        else
            count += 1
            i += 1
            j += 1
        end
    end
    return count
end

function jacquard_index(hash1, hash2)
    return hashcount(hash1, hash2) / min(length(hash1), length(hash2))
end

function mash_distance(hash1, hash2, k)
    j = jacquard_index(hash1, hash2)
    return -log(2 * j / (1 + j)) / k
end

function mash_distance(j, k)
    return -log(2 * j / (1 + j)) / k
end
