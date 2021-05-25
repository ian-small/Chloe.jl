include("Utilities.jl")

function embl2sff(emblfilepath::String, sfffilepath::String, seqfilepath::String)
    open(emblfilepath, "r") do f
        outfile = open(sfffilepath, "w")

        #ID and length
        idline = readline(f)
        @assert startswith(idline, "ID   ")
        tags = split(idline[6:end], ";")
        ID = tags[1]
        lengthtags = split(last(tags))
        seqlength = parse(Int32, lengthtags[1])

        seqfile = open(seqfilepath, "w")
        write(seqfile, ">", ID, "\n")

        #Features
        fmodels = Vector{Vector{Feature}}()
        rmodels = Vector{Vector{Feature}}()

        local have_model::Bool, waiting_for_coords::Bool, gene::String, product::String, type::String, coords::String
        function reset()
            have_model = false
            waiting_for_coords = false
            gene = ""
            product = ""
            type = ""
            coords = ""
        end

        reset()
        readingseq = false
        while !eof(f)
            line = readline(f)
            if startswith(line, "FT")
                maybetype = split(line[6:21])
                if !isempty(maybetype)
                    if have_model && type in ["CDS", "tRNA", "rRNA"]
                        #println(gene, '\t', product)
                        if haskey(gene_names, gene); gene = gene_names[gene];
                        elseif haskey(gene_names, product); gene = gene_names[product]; end
                        if gene == ""; gene = product; end
                        #println(gene, '\t', product)
                        (fmodel, rmodel) = embl2model(gene, type, coords, seqlength)
                        !isempty(fmodel) && append!(fmodels, cleanupfeatures(fmodel, seqlength))
                        !isempty(rmodel) && append!(rmodels, cleanupfeatures(rmodel, seqlength))
                        reset()
                    end
                    type = maybetype[1]
                    coords = line[22:end]
                    waiting_for_coords = true
                elseif startswith(line[22:end], "/")
                    waiting_for_coords = false
                    if startswith(line[22:end], "/gene=")
                        gene = replace(filter(x -> !isspace(x), line[29:end-1]), "_"=>"-")
                        have_model = true
                    elseif startswith(line[22:end], "/product=")
                        product = replace(filter(x -> !isspace(x), line[32:end-1]), "_"=>"-")
                        have_model = true
                    end
                elseif waiting_for_coords && (startswith(line[22:end], r"[0-9]") || startswith(line[22:end], "complement"))
                    coords *= line[22:end]
                end
            elseif startswith(line, "SQ")
                readingseq = true
            elseif startswith(line, "//")
                readingseq = false
            elseif readingseq
                write(seqfile, uppercase(join(split(line)[1:end-1])))
            end
        end
        println("fmodels: $(length(fmodels)) rmodels: $(length(rmodels))")
        addintrons!(fmodels)
        addintrons!(rmodels)
        writeSFF(outfile, ID, seqlength, fmodels, rmodels)
        close(outfile)
    end

end

mutable struct Feature
    gene::String
    type::String
    order::UInt8
    start::Int32
    length::Int32
    # phase is the number of nucleotides to skip at the start of the sequence 
    # to be in the correct reading frame
    phase::Int8
    function Feature(path, s, l, p)
        tags = split(path, "/")
        t = length(tags) == 4 ? tags[3] : tags[2]
        o = length(tags) == 4 ? parse(UInt8, tags[4][1]) : parse(UInt8, tags[3][1])
        return new(tags[1], t, o, s, l, p)
    end
end

const annotation_path(f::Feature) = begin
    join([f.gene, f.type, string(f.order)], "/")
end

#check if features overlap the end of the genome and should be merged
#split rps12 models into A and B models
function cleanupfeatures(model::Vector{Feature}, genome_length::Int32)
    models = Vector{Vector{Feature}}()
    firstf = first(model)
    lastf = last(model)
    if lastf.start + lastf.length == genome_length + 1
        if lastf.gene == firstf.gene
            if firstf.start == 1
            #must be contiguous and same feature if we got as far as here
            lastf.length += firstf.length
            deleteat!(model, 1)
            end
        end
    end
    #remove stop codon from CDS features
    if lastf.type == "CDS"
        lastf.length -= 3
    end
    if firstf.gene == "rps12A" && lastf.gene == "rps12B"
        rps12A = popfirst!(model)
        push!(models, [rps12A])
    end
    push!(models, model)
    return models
end

function encloses(text::String, start_index::Integer)
    enclosures = 0
    for i in start_index:length(text)
        if text[i] == '('
            enclosures += 1
        elseif text[i] ==')'
            enclosures -= 1
        end
        if enclosures == 0
            return i < length(text) ? false : true
        end
    end
    @assert false "broken enclosure: $text"
end

function embl2model(gene::AbstractString, type::AbstractString, coords::AbstractString, seqlength::Int32, strand::Char = '+')
    #println(coords)
    ffeatures = Vector{Feature}()
    rfeatures = Vector{Feature}()
    if startswith(coords, "complement") && encloses(coords, 11)
        ffeatures, rfeatures = embl2model(gene, type, coords[12:end-1], seqlength, '-')
    elseif startswith(coords, "join") && encloses(coords, 5)
        ffeatures, rfeatures = embl2model(gene, type, coords[6:end-1], seqlength, strand)
    elseif startswith(coords, "order") && endswith(coords, ")")
        ffeatures, rfeatures = embl2model(gene, type, coords[7:end-1], seqlength, strand)
    else
        ranges = split(coords, ",")
        if strand == '-'
            reverse!(ranges)
        end
        cumulative_length = 0
        for (index, range) in enumerate(ranges)
            #println(range)
            if gene == "rps12" || gene == "rps12A"
                gene = index == 1 ? "rps12A" : "rps12B"
            end
            (feature, fstrand) = embl2feature(gene, type, index, range, cumulative_length, seqlength, strand)
            if fstrand == '+'; push!(ffeatures, feature)
            elseif fstrand == '-'; push!(rfeatures, feature); end
            cumulative_length += feature.length
        end
    end
    return ffeatures, rfeatures
end

function embl2feature(gene::AbstractString, type::AbstractString, index::Int, range::AbstractString, cumulative_length::Int, seqlength::Int32, strand::Char)
    if startswith(range, "complement")
        strand = '-'
        range = range[12:end-1]
    end
    positions = split(range, r"[<>.]+", keepempty=false)
    start::Int32 = parse(Int32, positions[1])
    if length(positions) == 2
        stop = parse(Int32, positions[2])
    else
        stop = start
    end
    feature_length::Int32 =  stop - start + 1
    #reverse complement if - strand
    start = strand == '-' ? seqlength - stop + 1 : start
    phase = type == "CDS" ? phase_counter(Int8(0), cumulative_length) : 0
    order = 2 * index - 1
    if gene == "rps12B"; order -= 1; end;
    fname = gene*"/"*type*"/"*string(order)
    return Feature(fname, start, feature_length, phase), strand
end

function writeSFF(outfile::IO, id::AbstractString, genome_length::Int32, f_models::Vector{Vector{Feature}}, r_models::Vector{Vector{Feature}})
    
    function getmodelID!(model_ids::Dict{String,Int32}, model::Vector{Feature})
        gene_name = model[1].gene
        instance_count::Int32 = 1
        model_id = "$(gene_name)/$(instance_count)"
        while get(model_ids, model_id, 0) ≠ 0
            instance_count += 1
            model_id = "$(gene_name)/$(instance_count)"
        end
        model_ids[model_id] = instance_count
        return model_id
    end

    function out(outfile::IO)
        model_ids = Dict{String,Int32}()
        write(outfile, id, "\t", string(genome_length), "\n")
        for model in f_models
            #println(model)
            isempty(model) && continue
            model_id = getmodelID!(model_ids, model)
            #println(model_id)
            write_model2SFF(outfile, model_id, model, '+')
        end
        for model in r_models
            isempty(model) && continue
            model_id = getmodelID!(model_ids, model)
            write_model2SFF(outfile, model_id, model, '-')
        end
    end

    out(outfile)
end

function write_model2SFF(outfile::IO, model_id::String, model::Vector{Feature}, strand::Char)

    for feature in model
        #println(feature)
        relative_length = get_relative_length(feature)
        #println(relative_length)
        if relative_length == 0
            continue #avoid writing unrecognised features
        end
        write(outfile, model_id)
        write(outfile, "/")
        write(outfile, feature.type)
        write(outfile, "/")
        write(outfile, string(feature.order))
        write(outfile, "\t")
        write(outfile, join([strand,string(feature.start),string(feature.length),string(feature.phase)], "\t"))
        write(outfile, "\t", @sprintf("%.3f",relative_length))
        write(outfile, "\t0\t\t") #0 for stack depth, empty notes column
        write(outfile, "\n")
    end

end

const unknown_features = Dict{String, Int}()

function get_relative_length(f::Feature)
    fname = annotation_path(f)
    template = get(templates, fname, nothing)
    if isnothing(template)
        #println("failed to find template $fname")
        count = 0
        if haskey(unknown_features, fname); count = unknown_features[fname]; end
        unknown_features[fname] = count + 1
        return 0
    end
    return f.length/template.median_length
end

struct FeatureTemplate
    path::String  # similar to .sff path
    median_length::Float32
end

const templates = Dict{String,FeatureTemplate}()

function read_templates(file::String)
    if !isfile(file)
        error("\"$(file)\" is not a file")
    end
    if filesize(file) === 0
        error("no data in \"$(file)!\"")
    end
    open(file) do f
        readline(f) #skip header
        for line in eachline(f)
            fields = split(line, '\t')
            template = FeatureTemplate(fields[1], parse(Float32, fields[3]))
            templates[template.path] = template
        end
    end
end

const gene_names = Dict{String, String}()

function read_name_mappings(file::String)
    open(file) do f
        for line in eachline(f)
            names = split(line,"\t")
            gene_names[names[1]] = names[2]
        end
    end
end

function addintrons!(models::Vector{Vector{Feature}})
    for model in models
        exons = length(model)
        i = 1
        while i < exons
            #insert intron
            intron_order = i + 1
            if model[1].gene == "rps12B"; intron_order += 1; end
            intronpath = join([model[1].gene, "intron", string(intron_order)], "/")
            preceding_exon = model[i * 2 - 1]
            insert!(model, i * 2, Feature(intronpath, preceding_exon.start + preceding_exon.length, model[i * 2].start - preceding_exon.start - preceding_exon.length, 0))
            i += 1
        end
    end
end

function dir2sff(dir::String)
    files = filter(x->endswith(x, ".embl"), readdir(dir, join = true))
	for file in files
		id = split(basename(file),".")[1]
        println(id)
        outfile = joinpath(dir, id*".sff")
        seqfile = joinpath(dir, id*".fasta")
		embl2sff(file, outfile, seqfile)
	end
    for (f, c) in sort(collect(unknown_features), by = x -> x[2], rev = true)
        println(f, "\t", c)
    end
end