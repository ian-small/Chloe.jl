module Sff2Gff
export writeallGFF3

import ..Annotator: Feature, SFF_Feature, readFeatures, groupFeaturesIntoGeneModels

struct ModelArray
    genome_id::String
    genome_length::Int32
    strand::Char
    features::Vector{SFF_Feature}
end
function mergeAdjacentFeaturesinModel!(model::Vector{SFF_Feature}, genome_id, genome_length, strand)
    f1_index = 1
    f2_index = 2
    while f2_index <= length(model)
        f1 = model[f1_index].feature
        f2 = model[f2_index].feature
        # if adjacent features are same type, merge them into a single feature
        if f1.type == f2.type && f2.start - f1.start - f1.length ≤ 100
            @debug "[$(genome_id)]$(strand) merging adjacent $(f1.path) and $(f2.path)"
            f1.length += f2.length
            deleteat!(model, f2_index)
        else
            f1_index += 1
            f2_index += 1
        end
    end
    return ModelArray(genome_id, genome_length, strand, model)
end

function sff2gffcoords(f::Feature, strand::Char, genome_length::Integer)
    start = f.start
    finish = f.start + f.length - 1
    length = finish - start + 1
    
    if strand == '-'
        start = genome_length - finish + 1
        if start < 1
            start += genome_length
        end
        finish = start + length - 1
    end
    return (start, finish, length)
end

function writeGFF3(outfile, genemodel::ModelArray; sep::Bool=true)
    
    function write_line(type, start, finish, id, parent, phase="."; key="Parent")
        l = [genemodel.genome_id, "Chloe", type, start, finish, ".", genemodel.strand, phase]
        write(outfile, join(l, "\t"))
        write(outfile, "\t", "ID=", id, ";", key, "=", parent, "\n")
    end

    features = genemodel.features
    id = first(features).feature.gene
    id = path_components[1]
    if parse(Int, path_components[2]) > 1
        id = id * "-" * path_components[2]
    end

    start = minimum(f.start for f in features)

    if ft == "CDS" && !startswith(id, "rps12A")
        last(features).length += 3 #add stop codon
    end
    finish = maximum(f.start + f.length - 1 for f in features)
    ft = getFeatureType(first(features))

    length = finish - start + 1
    
    if genemodel.strand == '-'
        start = genemodel.genome_length - finish + 1
        if start < 1
            start += genemodel.genome_length
        end
        finish = start + length - 1
    end

    # gene
    write_line("gene", start, finish, id, path_components[1]; key="Name")
    # RNA product
    parent = id
    if ft == "CDS"
        parent =  id * ".mRNA"
        write_line("mRNA", start, finish, parent, id)
    elseif ft == "rRNA"
        parent = id * ".rRNA"
        write_line("rRNA", start, finish, parent, id)
    elseif ft == "tRNA"
        parent =  id * ".tRNA"
        write_line("tRNA", start, finish, parent, id)
    end
    for feature in (genemodel.strand == '+' ? features : Iterators.reverse(features))
        type = getFeatureType(feature)
        if type == "tRNA" || type == "rRNA"
            type = "exon"
        end
        
        start, finish, length = wrap(feature, genemodel.strand, genemodel.genome_length)
        
        phase = type == "CDS" ? string(feature.phase) : "."
        write_line(type, start, finish, feature.path, parent, phase)
    end
    if sep
        write(outfile, "###\n")
    end
end

function writeallGFF3(;sff_files=String[], directory=nothing, sep::Bool=true)
    
    function add_models!(models_as_feature_arrays, features, strand)
        # afeatures = collect(features.interval_tree, Vector{Feature}())
        afeatures =  Vector{SFF_Feature}()
        for feature in features.interval_tree
            push!(afeatures, SFF_Feature(feature, 0.0, 0.0, 0.0, 0.0, 0.0))
        end
        sort!(afeatures, by=f -> [f.feature.gene, f.feature.order])
        models = groupFeaturesIntoGeneModels(afeatures)
        for model in models
            if isempty(model)
                continue
            end
            push!(models_as_feature_arrays, mergeAdjacentFeaturesinModel!(model, features.genome_id, features.genome_length, strand))
        end
    end
    
    for infile in sff_files
        features = readFeatures(infile)
        models_as_feature_arrays = Vector{ModelArray}()
        # for each strand group into gene models
        add_models!(models_as_feature_arrays, features.forward, '+')
        add_models!(models_as_feature_arrays, features.reverse, '-')

        # interleave gene models from both strands
        sort!(models_as_feature_arrays, by=m -> m.features[1].feature.start)

        # write models in GFF3 format
        fname = features.forward.genome_id * ".gff3";
        if directory !== nothing
            fname = joinpath(directory, fname)
        else
            d = splitpath(infile)
            fname = joinpath(d[1:end - 1]..., fname)
        end
        @info "writing gff3: $fname"
        open(fname, "w") do outfile
            write(outfile, "##gff-version 3.2.1\n")
            for model in models_as_feature_arrays
                writeGFF3(outfile, model; sep=sep)
            end
        end
    end
end
end # module
