
using BioSequences
using GenomicAnnotations
using UUIDs

function write_model2SFF(outfile::IO, model::SFF_Model)
    isnothing(model) && return
    model_id = model.gene * "/" * string(model.gene_count)
    for sff in model.features
        f = sff.feature
        id = "$(model.gene)/$(model.gene_count)/$(f.type)/$(f.order)"
        write(outfile, id)
        write(outfile, "\t")
        write(outfile, join([model.strand, string(f.start), string(f.length), string(f.phase)], "\t"))
        write(outfile, "\t")
        write(
            outfile,
            join(
                [
                    @sprintf("%.3g", sff.relative_length),
                    @sprintf("%.3g", sff.stackdepth),
                    @sprintf("%.3g", sff.gmatch),
                    @sprintf("%.3g", sff.feature_prob),
                    @sprintf("%.3g", sff.coding_prob)
                ],
                "\t"
            )
        )
        write(outfile, "\t")
        write(outfile, join(model.warnings, "; "))
        write(outfile, "\n")
    end
    for warning in model.warnings
        @warn("$(model_id) $(warning)")
    end
end

function writeSFF(
    outfile::Union{String,IO},
    id::String, # NCBI id
    genome_length::Int32,
    mean_coverage::Float32,
    models::FwdRev{Vector{SFF_Model}}
)
    function out(outfile::IO)
        write(outfile, id, "\t", string(genome_length), "\t", @sprintf("%.3f", mean_coverage), "\n")
        for model in models.forward
            write_model2SFF(outfile, model)
        end
        for model in models.reverse
            write_model2SFF(outfile, model)
        end
    end
    if outfile isa String
        maybe_gzwrite(outfile::String) do io
            out(io)
        end
    else
        out(outfile)
    end
end

function construct_locus(sff::SFF_Model, feature_type::String, target_length::Integer)
    sfeatures = filter(x -> x.feature.type == feature_type, sff.features)
    isempty(sfeatures) && return nothing
    loci = Vector{AbstractLocus}()
    for sfeature in sfeatures
        f = sfeature.feature
        start = sff.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, target_length)
        span = ClosedSpan(start:start + f.length-1)
        push!(loci, sff.strand == '+' ? span : Complement(span))
    end
    if sff.strand == '-'; reverse!(loci); end
    locus = length(loci) == 1 ? first(loci) : Join(loci)
    locus
end

function chloe2biojulia(chloe::ChloeAnnotation)::GenomicAnnotations.Record
    biojulia = GenomicAnnotations.Record()
    biojulia.name = chloe.target_id
    biojulia.header = "##gff-version 3\n##source-version Chloe 1.0\n##sequence-region\tchloe.target_id\t1\t$(chloe.target_length)\n"
    biojulia.circular = true
    # add

    sffs = vcat(chloe.annotation.forward, chloe.annotation.reverse)
    for sff in sffs
        merge_adjacent_features!(sff)
        startswith(sff.gene, "rps12") && continue
        ft = featuretype(sff)
        if ft ≠ "repeat_region"; ft = "gene"; end
        # construct gene/repeat_region feature
        span = gene_span(sff)
        if sff.strand == '-'; span = reverse_complement(span, chloe.target_length); end
        locus = ClosedSpan(span)
        if sff.strand == '-'; locus = Complement(locus); end
        gene_id = string(uuid4())
        addgene!(biojulia, Symbol(ft), locus; locus_tag = gene_id, ID = gene_id, gene = sff.gene, name = sff.gene)
        # optionally construct mRNA feature
        # construct CDS, tRNA or rRNA feature
        for feature_type in ["CDS", "tRNA", "rRNA"]
            locus = construct_locus(sff, feature_type, chloe.target_length)
            if ~isnothing(locus)
                id = string(uuid4())
                addgene!(biojulia, Symbol(feature_type), locus; parent = gene_id, locus_tag = id, ID = id, name = "$(sff.gene).$feature_type")
            end
        end
        # construct intron feature(s)
        introns = filter(x -> x.feature.type == "intron", sff.features)
        for (i, intron) in enumerate(introns)
            f = intron.feature
            start = sff.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
            locus = ClosedSpan(start:start + f.length-1)
            if sff.strand == '-'; locus = Complement(locus); end
            id = string(uuid4())
            addgene!(biojulia, :intron, locus; parent = gene_id, locus_tag = id, ID = id, name = "$(sff.gene).intron.$i" )
        end
    end
    # join rps12A and rps12B features
    for a in filter(x -> x.gene == "rps12A", sffs), b in filter(x -> x.gene == "rps12B", sffs)
        # construct gene/repeat_region feature
        aspan = gene_span(a)
        if a.strand == '-'; aspan = reverse_complement(aspan, chloe.target_length); end
        alocus = ClosedSpan(aspan)
        if a.strand == '-'; alocus = Complement(alocus); end
        bspan = gene_span(b)
        if b.strand == '-'; bspan = reverse_complement(bspan, chloe.target_length); end
        blocus = ClosedSpan(bspan)
        if b.strand == '-'; blocus = Complement(blocus); end
        gene_id = string(uuid4())
        addgene!(biojulia, :gene, Join([alocus, blocus]); locus_tag = gene_id, ID = gene_id, gene = "rps12", name = "rps12")
        # construct CDS feature
        alocus = construct_locus(a, "CDS", chloe.target_length)
        blocus = construct_locus(b, "CDS", chloe.target_length)
        id = string(uuid4())
        addgene!(biojulia, :CDS, Join([alocus, blocus]); parent = gene_id, locus_tag = id, ID = id, name = "rps12.CDS")
        # construct rps12A internal intron feature(s) (I don't think there are any, but just in case...)
        aintrons = filter(x -> x.feature.type == "intron", a.features)
        intron_count = 1
        for intron in aintrons[1:end-1]
            f = intron.feature
            start = a.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
            locus = ClosedSpan(start:start + f.length-1)
            if a.strand == '-'; locus = Complement(locus); end
            id = string(uuid4())
            addgene!(biojulia, :intron, locus; parent = gene_id, locus_tag = id, ID = id, name = "rps12.intron.$intron_count" )
            intron_count += 1
        end
        bintrons = filter(x -> x.feature.type == "intron", b.features)
        # construct transpliced intron feature
        f = last(aintrons).feature
        start = a.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
        alocus = ClosedSpan(start:start + f.length-1)
        if a.strand == '-'; alocus = Complement(alocus); end
        f = first(bintrons).feature
        start = b.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
        blocus = ClosedSpan(start:start + f.length-1)
        if b.strand == '-'; blocus = Complement(blocus); end
        id = string(uuid4())
        addgene!(biojulia, :intron, Join([alocus, blocus]); parent = gene_id, locus_tag = id, ID = id, name = "rps12.intron.$intron_count" )
        intron_count += 1
        # construct rps12B internal intron feature(s)
        for intron in bintrons[2:end]
            f = intron.feature
            start = b.strand == '+' ? f.start : reverse_complement(f.start + f.length-1, chloe.target_length)
            locus = ClosedSpan(start:start + f.length-1)
            if b.strand == '-'; locus = Complement(locus); end
            id = string(uuid4())
            addgene!(biojulia, :intron, locus; parent = gene_id, locus_tag = id, ID = id, name = "rps12.intron.$intron_count" )
            intron_count += 1
        end
    end
    sort!(biojulia.genes)
    biojulia
end

function write_result(config::ChloeConfig, target::FwdRev{CircularSequence}, result::ChloeAnnotation, filestem::String)::Tuple{Union{String,IO},String}
    if ~config.no_transform
        FASTAWriter(open(filestem * ".chloe.fa", "w")) do outfile
            write(outfile, FASTARecord(result.target_id, target.forward[1:length(target.forward)]))
        end
    end
    if config.sff
        out = filestem * ".chloe.sff"
        writeSFF(out, result.target_id, result.target_length, geomean(values(result.coverages)), result.annotation)
    end
    if ~config.no_gff || config.gbk || config.embl
        biojulia = chloe2biojulia(result)
        #linearise
        biojulia.sequence = target.forward[1:length(target.forward)]
        if ~config.no_gff
            out = filestem * ".chloe.gff"
            GFF.printgff(out, biojulia)
        end
        if config.gbk
            out = filestem * ".chloe.gbk"
            GenBank.printgbk(out, biojulia)
        end
        if config.embl
            out = filestem * ".chloe.embl"
            EMBL.printembl(out, biojulia)
        end
    end
    return out, result.target_id
end
