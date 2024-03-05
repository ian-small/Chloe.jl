const startcodon = biore"ATG"d
const stopcodon = biore"(TAG)|(TAA)|(TGA)"d

struct FramedMatches
    genome_length::Int32
    frames::Vector{CircularVector}
end

function codonmatches(seq::BioSequences.NucleicSeqOrView, pattern)
    frames = Vector{CircularVector}(undef, 3)
    for f in 1:3
        frames[f] = CircularVector(Vector{Int32}())
    end
    for m in eachmatch(pattern, seq)
        i::Int32 = m.captured[1]
        push!(frames[mod1(i, 3)], i)
    end
    return FramedMatches(length(seq), frames)
end

function framedstops(seq::CircularSequence)
    matchmap = codonmatches(seq[1:seq.length], stopcodon)
    penultimate::Int32 = seq.length - one(Int32)
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
    matchmap = codonmatches(seq[1:seq.length], startcodon)
    penultimate::Int32 = seq.length - one(Int32)
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
    frame = mod1(i, 3)
    framestops = orfmap.stops[frame]
    #find first stop after index and last stop before
    index::Int32 = searchsortedfirst(framestops.v, i)
    stop = framestops[index]
    pstop = framestops[index-one(Int32)]
    if stop > pstop
        return pstop, stop - pstop + 3
    else
        return pstop, orfmap.genome_length - pstop + stop + 3
    end
end

function getmaxorflength(orfmap::FramedMatches, index)
    maxorf = 0
    for i in index-1:index+1
        orf = getorfcoords(orfmap, i)
        if orf[2] > maxorf
            maxorf = orf[2]
        end
    end
    return maxorf
end

function isstartcodon(codon::Union{Nothing,Tuple{DNA,DNA,DNA}}, allow_editing::Bool=false, allow_GTG::Bool=false)::Bool
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

function isstopcodon(codon::Union{Nothing,Tuple{DNA,DNA,DNA}}, allow_editing::Bool=false)::Bool
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
    return false
end

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
                push!(orfs, Feature("unassigned_orf_" * strandprefix * string(orf_counter) * "/CDS/1", Int32(pstop + 3), Int32(distance_from_pstop - 3), Int8(0)))
                orf_counter += 1
            end
            pstop = stop
        end
    end
    # add orfs that cross end of genome
    for f1 in stops.frames, f2 in stops.frames
        dist = circulardistance(last(f1.v), first(f2.v), length(seq))
        if dist % 3 == 0 #in frame
            if dist ≥ minorflength
                push!(orfs, Feature("unassigned_orf_" * strandprefix * string(orf_counter) * "/CDS/1", Int32(last(f1.v) + 3), Int32(dist - 3), Int8(0)))
            end
            continue
        end
    end
    return sort!(orfs, by=x -> x.length, rev=true)
end

const exon_indices = Dict{Int,String}()
open(joinpath(@__DIR__, "exon_indices.txt")) do exons
    for line in readlines(exons)
        tags = split(line, "\t")
        exon_indices[parse(Int, tags[1])] = tags[2]
    end
end

Codon = Tuple{DNA,DNA,DNA}

const codons = Dict{Codon,Int16}()
for i in [DNA_A, DNA_C, DNA_G, DNA_T], j in [DNA_A, DNA_C, DNA_G, DNA_T], k in [DNA_A, DNA_C, DNA_G, DNA_T]
    codons[(i, j, k)] = 0
end
const codon_strings = Vector{String}(undef, length(codons))
for (i, k) in enumerate(sort(collect(keys(codons))))
    codon_strings[i] = String([convert(Char, k[1]), convert(Char, k[2]), convert(Char, k[3])])
end

function countcodons(orf::Feature, seq::CircularSequence)::Vector{Float64}
    #initialise codon counts
    for (k, v) in codons
        codons[k] = 0
    end
    stopcount = orf.start + orf.length - 2
    codonstart::Int32 = orf.start + orf.phase
    while codonstart < stopcount
        codon = (seq[codonstart], seq[Int32(codonstart + 1)], seq[Int32(codonstart + 2)])
        if haskey(codons, codon)
            count = codons[codon]
            codons[codon] = count + 1
        end
        codonstart += 3
    end
    #convert counts to frequencies
    sortedcounts = sort(collect(codons))
    countsum = sum(map(x -> x[2], sortedcounts))
    codonfrequencies::Vector{Float64} = map(x -> x[2] / countsum, sortedcounts)
    # add the relative length of the ORF; assuming no ORFs are > 3000 codons this should be between 0-1
    push!(codonfrequencies, countsum / 3000)
    return codonfrequencies
end
const xgb_coding_cache = XGBoost.DMatrix[]
const xgb_coding_model = XGBoost.Booster(xgb_coding_cache, model_file=joinpath(@__DIR__, "xgb.coding.model"))
function xgb_coding_classifier(codonfrequencies::Vector{Float64})::Float32
    #explicit test for stop codons
    codonfrequencies[49] > 0 && return Float32(0.0) #TAA
    codonfrequencies[51] > 0 && return Float32(0.0) #TAG
    codonfrequencies[57] > 0 && return Float32(0.0) #TGA
    #test with XGB (is there a better way to construct the input matrix?)
    pred = XGBoost.predict(xgb_coding_model, reshape(codonfrequencies[vcat(collect(1:48), [50], collect(52:56), collect(58:end))], 1, :))
    return pred[1]
end
