const startcodon = biore"ATG"d
const stopcodon = biore"(TAG)|(TAA)|(TGA)"d

struct FramedMatches
    genome_length::Int32
    frames::Vector{CircularVector}
end

function codonmatches(seq::LongDNASeq, pattern)
    frames = Vector{CircularVector}(undef,3)
    for f in 1:3
        frames[f] = CircularVector(Vector{Int32}())
    end
    for m in eachmatch(pattern, seq)
        i::Int32 = m.captured[1]
        push!(frames[mod1(i,3)],i)
    end
    return FramedMatches(length(seq), frames)
end

function framedstops(seq::CircularSequence)
    matchmap = codonmatches(seq.sequence, stopcodon)
    penultimate::Int32 = length(seq.sequence) - one(Int32)
    penultimate_codon = getcodon(seq, penultimate)
    ultimate::Int32 = penultimate + one(Int32)
    ultimate_codon = getcodon(seq, ultimate)
    if isstopcodon(penultimate_codon)
        push!(matchmap.frames[mod1(penultimate, 3)], penultimate)
    end
    if isstopcodon(ultimate_codon)
        push!(matchmap.frames[mod1(ultimate, 3)], ultimate)
    end
    return matchmap
end

function framedstarts(seq::CircularSequence)
    matchmap = codonmatches(seq.sequence, startcodon)
    penultimate::Int32 = length(seq.sequence) - one(Int32)
    penultimate_codon = getcodon(seq, penultimate)
    ultimate::Int32 = penultimate + one(Int32)
    ultimate_codon = getcodon(seq, ultimate)
    if isstartcodon(penultimate_codon)
        push!(matchmap.frames[mod1(penultimate, 3)], penultimate)
    end
    if isstartcodon(ultimate_codon)
        push!(matchmap.frames[mod1(ultimate, 3)], ultimate)
    end
    return matchmap
end

function getorfcoords(orfmap::FramedMatches, i)
    frame = mod1(i,3)
    framestops = orfmap.stops[frame]
    #find first stop after index and last stop before
    index::Int32 = searchsortedfirst(framestops.v, i)
    stop = framestops[index]
    pstop = framestops[index - one(Int32)]
    if stop > pstop
        return pstop, stop-pstop + 3
    else
        return pstop, orfmap.genome_length - pstop + stop + 3
    end
end

function getmaxorflength(orfmap::FramedMatches, index)
    maxorf = 0;
    for i in index-1:index+1
        orf = getorfcoords(orfmap, i)
        if orf[2] > maxorf; maxorf = orf[2]; end
    end
    return maxorf
end

function isstartcodon(codon::Union{Nothing,Tuple{DNA,DNA,DNA}}, allow_editing::Bool = false, allow_GTG::Bool = false)::Bool
    isnothing(codon) && return false
    if codon[1] == DNA_A || (allow_GTG && codon[1] == DNA_G)
        if codon[2] == DNA_T || (allow_editing && codon[2] == DNA_C)
            if codon[3] == DNA_G #ATG / GTG / ACG (GCG)
                return true
            end
        end
    end
    return false
end

function isstopcodon(codon::Union{Nothing,Tuple{DNA,DNA,DNA}}, allow_editing::Bool = false)::Bool
    isnothing(codon) && return false
    if codon[1] == DNA_T || (allow_editing && codon[1] == DNA_C)
        if codon[2] == DNA_A
            if codon[3] == DNA_A || codon[3] == DNA_G    #TAA, TAG / CAA, CAG
                return true
            end
        elseif codon[2] == DNA_G && codon[3] == DNA_A  #TGA / CGA
            return true
        end
    end
    #println("false")
    return false;
end

#= function getallorfs(seq::CircularSequence, strand::Char, minorflength::Int32)
    starts = framedstarts(seq)
    stops = framedstops(seq)
    orfs = AFeature()
    orf_counter = 1
    strandprefix = strand == '+' ? "f" : "r"
    for f in 1:3
        frame = stops.frames[f]
        pstop = frame[Int32(-1)]
        for stop in frame[1:length(frame)]
            distance_from_pstop = circulardistance(pstop, stop, length(seq))
            if distance_from_pstop ≥ minorflength
                #worth looking for start
                start = starts.frames[f][Int32(searchsortedfirst(starts.frames[f].v, pstop))]
                distance_from_start = circulardistance(start, stop, length(seq))
                if distance_from_start ≥ minorflength && distance_from_start ≤ distance_from_pstop
                    #println(pstop, '\t', start, '\t', stop, '\t', distance_from_start)
                    push!(orfs, Feature("unassigned_orf_" * strandprefix * string(orf_counter) * "/1/CDS/1", start, distance_from_start + 3, 0))
                    orf_counter += 1
                end
            end
            pstop = stop
        end
    end
    return orfs
end
 =#
function getallorfs(seq::CircularSequence, strand::Char, minorflength::Int32)
    stops = framedstops(seq)
    orfs = Feature[]
    orf_counter = 1
    strandprefix = strand == '+' ? "f" : "r"
    for f in 1:3
        frame = stops.frames[f]
        pstop = frame[Int32(1)]
        for stop in frame[2:length(frame)]
            distance_from_pstop = circulardistance(pstop, stop, length(seq))
            if distance_from_pstop ≥ minorflength
                push!(orfs, Feature("unassigned_orf_" * strandprefix * string(orf_counter) * "/CDS/1", pstop+3, distance_from_pstop-3, 0))
                orf_counter += 1
            end
            pstop = stop
        end
    end
    #add orfs that cross end of genome
    for f1 in stops.frames, f2 in stops.frames
        dist = circulardistance(last(f1.v), first(f2.v), length(seq))
        if dist % 3 == 0 #in frame
            if dist ≥ minorflength
                push!(orfs, Feature("unassigned_orf_" * strandprefix * string(orf_counter) * "/CDS/1", last(f1.v)+3, dist-3, 0))
            end
            continue
        end
    end
    return sort!(orfs, by = x -> x.length, rev = true)
end

const exon_indices = Dict{Int, String}()
open(joinpath(@__DIR__, "exon_indices.txt")) do exons
    for line in readlines(exons)
        tags = split(line, "\t")
        exon_indices[parse(Int, tags[1])] = tags[2]
    end
end

Codon = Tuple{DNA,DNA,DNA}

const codons = Dict{Codon, Int16}()
for i in [DNA_A, DNA_C, DNA_G, DNA_T], j in [DNA_A, DNA_C, DNA_G, DNA_T], k in [DNA_A, DNA_C, DNA_G, DNA_T]
    codons[(i, j, k)] = 0
end
const codon_strings = Vector{String}(undef, length(codons))
for (i,k) in enumerate(sort(collect(keys(codons))))
    codon_strings[i] = String([convert(Char, k[1]), convert(Char, k[2]), convert(Char, k[3])])
end

function countcodons(orf::Feature, seq::CircularSequence)::Vector{Float64}
    #initialise codon counts
    for (k,v) in codons
        codons[k] = 0
    end
    stopcount = orf.start + orf.length - 2
    codonstart::Int32 = orf.start + orf.phase
    while codonstart < stopcount
        codon = (seq[codonstart], seq[Int32(codonstart+1)], seq[Int32(codonstart+2)])
        if haskey(codons, codon)
            count = codons[codon]
            codons[codon] = count + 1
        end
        codonstart += 3
    end
    #convert counts to frequencies
    sortedcounts = sort(collect(codons))
    countsum = sum(map(x -> x[2], sortedcounts))
    codonfrequencies::Vector{Float64} = map(x -> x[2]/countsum, sortedcounts)
    #add the relative length of the ORF; assuming no ORFs are > 3000 codons this should be between 0-1
    push!(codonfrequencies, countsum/3000)
    return codonfrequencies
end

function readGLMcodingclassifer()
    glm_coding_coeffs = Dict{String, Float64}()
    open(joinpath(@__DIR__, "GLMCodingClassifier.tsv")) do coding_coeffs
        readline(coding_coeffs) # skip header
        for line in readlines(coding_coeffs)
            fields = split(line, "\t")
            glm_coding_coeffs[fields[1]] = parse(Float64,fields[2])
        end
    end
    return glm_coding_coeffs
end

const GLM_Coding_Coeffs = readGLMcodingclassifer()

function glm_coding_classifier(codonfrequencies::Vector{Float64})::Float32
    #explicit test for stop codons
    codonfrequencies[49] > 0 && return Float32(0.0) #TAA
    codonfrequencies[51] > 0 && return Float32(0.0) #TAG
    codonfrequencies[57] > 0 && return Float32(0.0) #TGA
    #test with GLM
    pred = 0.0
    for (i,v) in enumerate(codonfrequencies)
        if i < 65
            key = codon_strings[i]
        else
            key = "orf_length"
        end
        coef = get(GLM_Coding_Coeffs, key, 0.0)
        pred += coef * v
    end
    pred += get(GLM_Coding_Coeffs, "(Intercept)", 0.0)
    odds = exp(pred)
    return Float32(odds / (1.0 + odds))
end
