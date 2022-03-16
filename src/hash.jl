using BioSequences
using FASTX
using Statistics

function sequence_entropy(seq::LongDNA)
    comp = zeros(Int, 4)
    for i::Int32 in 1:length(seq)
        if seq[i] == DNA_A
            comp[1] += 1
        elseif seq[i] == DNA_C
            comp[2] += 1
        elseif seq[i] == DNA_G
            comp[3] += 1
        elseif seq[i] == DNA_T
            comp[4] += 1
        end
    end
    entropy = 0.0
    for v in comp
        prob = v / length(seq)
        prob == 0 && continue
        entropy -= prob * log2(prob)
    end
    return comp, entropy
end

function minhash(seq::CircularSequence)
    numhashes = length(seq)
    hashes = UInt64[]
    sizehint!(hashes, numhashes)
    for i in 1:numhashes
        any(seq.mask[i:i+KMERSIZE-1]) && continue
        kmer = seq[i:i+KMERSIZE-1]
        push!(hashes, hash(canonical(kmer)))
    end
    return sort!(hashes)[1:SKETCHSIZE]
end

function minhash_references(; fasta_files=Vector{String}, output="reference_minhashes.hash")
    seqs = Vector{FASTA.Record}()
    for file in fasta_files
        println(file)
        if isdir(file)
            minhash_references(fasta_files=filter(x->endswith(x, ".fasta"), readdir(file, join = true)), output=output)
        end
        reader = open(FASTA.Reader, file)
        record = FASTA.Record() 
        read!(reader, record)
        push!(seqs, record)
        close(reader)
    end

    hashes = Dict{String,Vector{Int64}}()
    for seqrec in seqs
        seq = FASTA.sequence(LongDNA{4}, seqrec)
        hashes[FASTA.identifier(seqrec)] = minhash(CircularSequence(seq))
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
    #add sizehint! once number of references is established
    open(infile, "r") do hashfile
        while !eof(hashfile)
            length_id = read(hashfile, UInt8)
            id = String(read(hashfile, length_id, all=true))
            h = Vector{Int64}(undef,SKETCHSIZE)
            hashes[id] = read!(hashfile, h)
        end
    end
    return hashes
end

function searchhashes(hash, hashes)::Vector{Tuple{String,Int}}
    scores = Vector{Tuple{String,Int}}(undef,length(hashes))
    for (index, (id, h)) in enumerate(hashes)
        scores[index] = (id, hashcount(hash, h))
    end
    return sort!(scores, by = last, rev = true)
end

function hashcount(hash1, hash2)::Int
    l = length(hash1) + 1
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
    # hack to avoid picking a ref identical to the target
    return count == length(hash1) ? 0 : count
    #return count
end

function jacquard_index(hash1, hash2)
    return hashcount(hash1,hash2)/min(length(hash1),length(hash2))
end

function mash_distance(hash1, hash2, k)
    j = jacquard_index(hash1, hash2)
    return -log(2 * j / (1 + j))/k
end

function mash_distance(j, k)
    return -log(2 * j / (1 + j))/k
end
