include("Utilities.jl")

function embl2sff(filepath::String, outfilepath::String)
    open(filepath, "r") do f
        outfile = open(outfilepath, "w")

        #ID and length
        idline = readline(f)
        @assert startswith(idline, "ID   ")
        tags = split(idline[6:end], ";")
        ID = tags[1]
        lengthtags = split(last(tags))
        seqlength = parse(Int32, lengthtags[1])

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
        while !eof(f)
            line = readline(f)
            !startswith(line, "FT") && continue
            maybetype = split(line[6:21])
            if !isempty(maybetype)
                if have_model && type in ["CDS", "tRNA", "rRNA"]
                    if haskey(gene_names, gene); gene = gene_names[gene];
                    elseif haskey(gene_names, product); gene = gene_names[product]; end
                    if gene == ""; gene = product; end
                    (fmodel, rmodel) = embl2model(gene, product, type, coords, seqlength)
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
        end
        println("fmodels: $(length(fmodels)) rmodels: $(length(rmodels))")
        writeSFF(outfile, ID, seqlength, fmodels, rmodels)
        close(outfile)
    end

end

mutable struct Feature
    path::String
    start::Int32
    length::Int32
    # phase is the number of nucleotides to skip at the start of the sequence 
    # to be in the correct reading frame
    phase::Int8
    _path_components::Vector{String}
    Feature(path, start, length, phase) = new(path, start, length, phase, split(path, '/'))
end

function getFeatureName(feature::Feature)
    feature._path_components[1]
end

function annotationPath(feature::Feature)
    pc = feature._path_components
    join([pc[1],"?",pc[3], pc[4]], "/")
end

function getFeatureSuffix(feature::Feature)
    pc = feature._path_components
    join([pc[3], pc[4]], "/")
end

function rename(feature::Feature, new_name::String)
    feature.path = new_name
    feature._path_components = split(new_name, '/')
end

#check if features overlap the end of the genome and should be merged
#split rps12 models into A and B models
function cleanupfeatures(model::Vector{Feature}, genome_length::Int32)
    models = Vector{Vector{Feature}}()
    firstf = first(model)
    lastf = last(model)
    if lastf.start + lastf.length == genome_length + 1
        if getFeatureName(lastf) == getFeatureName(firstf)
            if firstf.start == 1
            #must be contiguous and same feature if we got as far as here
            lastf.length += firstf.length
            deleteat!(model, 1)
            end
        end
    end
    if getFeatureName(firstf) == "rps12A" && getFeatureName(lastf) == "rps12B"
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

function embl2model(gene::AbstractString, product::AbstractString, type::AbstractString, coords::AbstractString, seqlength::Int32, strand::Char = '+')
    #println(coords)
    ffeatures = Vector{Feature}()
    rfeatures = Vector{Feature}()
    if startswith(coords, "complement") && encloses(coords, 11)
        ffeatures, rfeatures = embl2model(gene, product, type, coords[12:end-1], seqlength, '-')
    elseif startswith(coords, "join") && encloses(coords, 5)
        ffeatures, rfeatures = embl2model(gene, product, type, coords[6:end-1], seqlength, strand)
    elseif startswith(coords, "order") && endswith(coords, ")")
        ffeatures, rfeatures = embl2model(gene, product, type, coords[7:end-1], seqlength, strand)
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
            (feature, fstrand) = embl2feature(gene, product, type, index, strand, range, cumulative_length, seqlength)
            if fstrand == '+'; push!(ffeatures, feature)
            elseif fstrand == '-'; push!(rfeatures, feature); end
            cumulative_length += feature.length
        end
    end
    return ffeatures, rfeatures
end

function embl2feature(gene::AbstractString, product::AbstractString, type::AbstractString, index::Int, strand::Char, range::AbstractString, cumulative_length::Int, seqlength::Int32)
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
    phase = type == "CDS" ? phaseCounter(Int8(0), cumulative_length) : 0
    if gene == "rps12B"; index -= 1; end;
    fname = gene*"/?/"*type*"/"*string(index)
    return Feature(fname, start, feature_length, phase), strand
end

function writeSFF(outfile::IO, id::AbstractString, genome_length::Int32, f_models::Vector{Vector{Feature}}, r_models::Vector{Vector{Feature}})
    
    function getModelID!(model_ids::Dict{String,Int32}, model::Vector{Feature})
        gene_name = getFeatureName(model[1])
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
            isempty(model) && continue
            model_id = getModelID!(model_ids, model)
            writeModelToSFF(outfile, model_id, model, '+')
        end
        for model in r_models
            isempty(model) && continue
            model_id = getModelID!(model_ids, model)
            writeModelToSFF(outfile, model_id, model, '-')
        end
    end

    out(outfile)
end

function writeModelToSFF(outfile::IO, model_id::String, model::Vector{Feature}, strand::Char)

    for feature in model
        relative_length = getRelativeLength(feature)
        relative_length == 0 && continue #avoid writing unrecognised features
        write(outfile, model_id)
        write(outfile, "/")
        write(outfile, getFeatureSuffix(feature))
        write(outfile, "\t")
        write(outfile, join([strand,string(feature.start),string(feature.length),string(feature.phase)], "\t"))
        write(outfile, "\t", @sprintf("%.3f",relative_length))
        write(outfile, "\t0\t\t") #0 for stack depth, empty notes column
        write(outfile, "\n")
    end

end

const unknown_features = Dict{String, Int}()

function getRelativeLength(f::Feature)
    fname = annotationPath(f)
    template = get(templates, fname, nothing)
    if isnothing(template)
        count = 0
        if haskey(unknown_features, fname); count = unknown_features[fname]; end
        unknown_features[fname] = count + 1
        return 0
    end
    return f.length/template.median_length
end

struct FeatureTemplate
    path::String  # similar to .sff path
    threshold_counts::Float32
    threshold_coverage::Float32
    median_length::Int32
end

const templates = Dict{String,FeatureTemplate}()

function readTemplates(file::String)::Dict{String,FeatureTemplate}

    if !isfile(file)
        error("\"$(file)\" is not a file")
    end
    if filesize(file) === 0
        error("no data in \"$(file)!\"")
    end

    gene_exons = String[]

    open(file) do f
        header = readline(f)
        for line in eachline(f)
            fields = split(line, '\t')
            path_components = split(fields[1], '/')
            if path_components[3] ≠ "intron"
                push!(gene_exons, path_components[1])
            end
            template = FeatureTemplate(fields[1], 
                parse(Float32, fields[2]),
                parse(Float32, fields[3]),
                parse(Int32, fields[4]))
            if haskey(templates, template.path)
                @error "duplicate path: $(template.path) in \"$file\""
            end
            templates[template.path] = template
            # push!(templates, template)
        end
    end
    return templates #, StatsBase.addcounts!(Dict{String,Int32}(), gene_exons)
end

const gene_names = Dict{String, String}()

function readNameMappings(file::String)
    open(file) do f
        for line in eachline(f)
            names = split(line,"\t")
            gene_names[names[1]] = names[2]
        end
    end
end

function dir2sff(dir::String)
    files = filter(x->endswith(x, ".embl"), readdir(dir, join = true))
	for file in files
		id = split(basename(file),".")[1]
        println(id)
        outfile = joinpath(dir, id*".sff")
        isfile(outfile) && continue
		embl2sff(file, outfile)
	end
    for (f, c) in unknown_features
        println(f, "\t", c)
    end
end