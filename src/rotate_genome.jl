using BioSequences
using FASTX

function rotategenome(infile::String, io::IO, flip_LSC::Bool, flip_SSC::Bool, extend::Int=0)
    records = open(infile) do io
        [record for record in FASTA.Reader(io)]
    end
    if isempty(records)
        error("unable to read sequence from $infile")

    elseif length(records) > 1
        @error "$infile contains multiple sequences; Chloë expects a single sequence per file"
    end
    fseq = CircularSequence(FASTA.sequence(LongDNA{4}, records[1]))
    rseq = reverse_complement(fseq)
    ir = Annotator.inverted_repeat(fseq, rseq)


    out = LongDNA(dna""d)

    if ir.blocklength >= 1000
        IR1 = (ir.src_index, ir.blocklength)
        IR2 = (genome_wrap(length(fseq), length(fseq) - ir.tgt_index - ir.blocklength + 2), ir.blocklength)
        scstart = genome_wrap(length(fseq), IR1[1] + ir.blocklength)
        SSC = (scstart, circulardistance(scstart, IR2[1], length(fseq)))
        scstart = genome_wrap(length(fseq), IR2[1] + ir.blocklength)
        LSC = (scstart, circulardistance(scstart, IR1[1], length(fseq)))
        if LSC[2] < SSC[2]
            temp = LSC
            LSC = SSC
            SSC = temp
            temp = IR1
            IR1 = IR2
            IR2 = temp
        end
        n = LSC[2] + IR1[2] + SSC[2] + IR2[2]
        @assert n == length(fseq) "$(FASTA.identifier(records[1])) $n != $(length(fseq))"
        lsc = LSC[1]:LSC[1]+LSC[2]-1
        ir1 = IR1[1]:IR1[1]+IR1[2]-1
        ssc = SSC[1]:SSC[1]+SSC[2]-1
        ir2 = IR2[1]:IR2[1]+IR2[2]-1
        @info "LSC=$(lsc)[$(length(lsc))] IR1=$(ir1)[$(length(ir1))] SSC=$(ssc)[$(length(ssc))] IR2=$(ir2)[$(length(ir2))]"
        append!(out, flip_LSC ? reverse_complement(fseq[lsc]) : fseq[lsc])
        append!(out, fseq[ir1])
        append!(out, flip_SSC ? reverse_complement(fseq[ssc]) : fseq[ssc])
        append!(out, fseq[ir2])
    else
        @warn "[$(infile)] no inverted repeat found: $(ir.blocklength) < 1000"
        append!(out, fseq)
    end

    if extend != 0
        e = min(extend > 0 ? extend : length(fseq) - 1, length(fseq) - 1)
        @debug "[$(infile)] extending by $e"
        # extend transformed seq
        append!(out, out[1:e])
    end

    writer = FASTA.Writer(io)
    id = FASTA.identifier(records[1]) * " rotated"
    rec = FASTA.Record(id, out)
    write(writer, rec)
    close(writer)
end

function rotategenomes(; fasta_files::Vector{String}, output::Union{String,Nothing}=nothing,
    flip_LSC::Bool=false, flip_SSC::Bool=false, extend::Int=0)

    for fasta in fasta_files
        d = splitpath(fasta)
        outfile = nothing
        if isnothing(output)
            io = stdout
        elseif isdir(output)
            f, ext = splitext(d[end])
            outfile = joinpath(output, f * "-rotated" * ext)
        else
            outfile = output
        end
        if !isnothing(outfile)
            io = open(outfile, "w")
            @info "$(fasta) rotating to $(outfile)"
        end
        try
            rotategenome(fasta, io, flip_LSC, flip_SSC, extend)
        finally
            if io != stdout
                close(io)
            end
        end
    end
end