const RNApairs = Dict((DNA_A, DNA_A) => 0, (DNA_A, DNA_C) => 0, (DNA_A, DNA_G) => 0, (DNA_A, DNA_T) => 1, (DNA_C, DNA_A) => 0, (DNA_C, DNA_C) => 0, (DNA_C, DNA_G) => 1, (DNA_C, DNA_T) => 0,
    (DNA_G, DNA_A) => 0, (DNA_G, DNA_C) => 1, (DNA_G, DNA_G) => 0, (DNA_G, DNA_T) => 1, (DNA_T, DNA_A) => 1, (DNA_T, DNA_C) => 0, (DNA_T, DNA_G) => 1, (DNA_T, DNA_T) => 0)


function complementaryRNAscore(seq1::LongDNASeq, seq2::LongDNASeq)
    @assert length(seq1) == length(seq2)
    score = 0
    for n in 1:length(seq1)
        score += get(RNApairs, (seq1[n], seq2[end - n + 1]), 0)
    end
    return score
end

function longest_complementary_match(seq1::LongDNASeq, seq2::LongDNASeq)
    rseq2 = reverse(seq2)
    len = length(rseq2)
    L = zeros(Int8, len, len)
    z = 0
    offsets = (0, 0)
    for i in 1:len, j in 1:len
        match = get(RNApairs, (seq1[i], rseq2[j]), 0)
        if match == 1
            #don't count terminal G:U pairs
            if (i == 1 || j == 1)
                if !((seq1[i] == DNA_G && rseq2[j] == DNA_T) || (seq1[i] == DNA_T && rseq2[j] == DNA_G))
                    L[i,j] = 1
                end
            elseif L[i - 1, j - 1] == 0
                if ((seq1[i] == DNA_G && rseq2[j] == DNA_T) || (seq1[i] == DNA_T && rseq2[j] == DNA_G))
                    L[i,j] = 0
                else
                    L[i,j] = 1
                end 
            else
                L[i, j] = L[i - 1, j - 1] + 1
            end
            if L[i, j] > z
                z = L[i, j]
                offsets = (i - z + 1, j - z + 1)
            end
        else
            L[i, j] = 0
        end
    end
    return offsets, z
end

function trimRNAends!(seq::CircularSequence, model::Vector{SFF_Feature})
    stem = 6
    slop = 3
    start = first(model).feature.start
    stop = last(model).feature.start + last(model).feature.length - 1
    offsets, stemlength = longest_complementary_match(seq[start - slop:start + slop + stem], seq[stop - stem - slop:stop + slop])
    if stemlength ≥ stem
        startchange = offsets[1] - slop - 1
        first(model).feature.start += startchange
        first(model).feature.length -= startchange
        stopchange = slop - offsets[2] + 1
        last(model).feature.length += stopchange
    end
end

function tRNAends!(seq::CircularSequence, model::Vector{SFF_Feature})
    start = first(model).feature.start
    stop = last(model).feature.start + last(model).feature.length - 1

    # before making any changes, check if model scores sufficiently highly to be left alone
    #could use template length, Tstemscore, acceptor stem score
    t = Tstemscore(seq[stop - 24 : stop])
    astart = first(model).feature.gene == "trnH-GUG" ? start + 1 : start
    a = complementaryRNAscore(seq[astart : astart + 6], seq[stop - 7 : stop - 1])
    println(first(model).feature.gene, "\t", t, "\t", a)
    t ≥ 14 &&  a ≥ 5 && return
    
    slop = 2
    maxscore = 0
    end3 = 0
    for i in stop - slop : stop + slop
        score = Tstemscore(seq[i - 24 : i])
        if score > maxscore
            maxscore = score
            end3 = i
        end
    end
    maxscore = 0
    end5 = 0
    for i in start - slop : start + slop
        score = complementaryRNAscore(seq[i : i + 6], seq[end3 - 7 : end3 - 1])
        if score > maxscore
            maxscore = score
            end5 = i
        end
    end
    if first(model).feature.gene == "trnH-GUG"; end5 -= 1; end #trnH-specific nucleotide at -1
    startchange = end5 - start
    first(model).feature.start += startchange
    first(model).feature.length -= startchange
    first(model).feature.start = genome_wrap(length(seq), first(model).feature.start)
    stopchange = stop - end3
    last(model).feature.length += stopchange

    t = Tstemscore(seq[stop - 24 : stop])
    astart = first(model).feature.gene == "trnH-GUG" ? start + 1 : start
    a = complementaryRNAscore(seq[astart : astart + 6], seq[stop - 7 : stop - 1])
    println(first(model).feature.gene, "\t", t, "\t", a)
end

function Tstemscore(seq::LongDNASeq)
    #relies on all chloroplast tRNAs having a canonical T-stem structure
    #scores 25 nt at 3' end of tRNA for match to TφC loop and T stem
    @assert length(seq) == 25
    score = 0
    #match to TφC
    if seq[6] == DNA_T && seq[7] == DNA_T && seq[8] == DNA_C
        score += 10
    end
    #T stem complementarity
    score += complementaryRNAscore(seq[1:5], seq[13:17])
    return score
end