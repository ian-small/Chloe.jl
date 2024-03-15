mutable struct Feature <: AbstractInterval{Int32}
    gene::String
    stack::CircularVector
    type::String
    order::UInt8
    start::Int32
    length::Int32
    # phase is the number of nucleotides to skip at the start of the sequence 
    # to be in the correct reading frame
    phase::Int8
    essential::Bool #from template
    median_length::Float32 #from template
end

function Feature(path::AbstractString, stack::CircularVector, s::Int32, l::Int32, p::Int8, e::Bool, ml::Float32)
    tags = split(path, "/")
    gene = tags[1]
    type = length(tags) == 4 ? tags[3] : tags[2]
    order = length(tags) == 4 ? parse(UInt8, tags[4][1]) : parse(UInt8, tags[3][1])
    return Feature(gene, stack, type, order, s, l, p, e, ml)
end

function Feature(path::AbstractString, s::Int32, l::Int32, p::Int8)
    tags = split(path, "/")
    gene = tags[1]
    type = length(tags) == 4 ? tags[3] : tags[2]
    order = length(tags) == 4 ? parse(UInt8, tags[4][1]) : parse(UInt8, tags[3][1])
    return Feature(gene, CircularVector(Int32[]), type, order, s, l, p, true, Float32(0))
end

Base.first(f::Feature) = f.start
function Base.last(f::Feature)
    f.start + f.length - one(Int32)
end

function Base.isless(a::Feature, b::Feature)
    return a.gene == b.gene ? isless(a.start + a.length / 2, b.start + b.length / 2) : isless(a.gene, b.gene)
end

const annotation_path(f::Feature) = begin
    join([f.gene, f.type, string(f.order)], "/")
end

AFeature = Vector{Feature}
AAFeature = Vector{AFeature}

struct FeatureTree
    ftree::IntervalTree{Int32,Feature}
    length::Int32
    wrapped_intervals::Dict{Feature,Feature}
end

function Base.push!(ctree::FeatureTree, interval::Feature)
    if last(interval) ≤ ctree.length
        push!(ctree.ftree, interval)
    else
        first_interval = Feature(interval.gene, interval.stack, interval.type, interval.order, first(interval), ctree.length - first(interval) + one(Int32), zero(Int8), interval.essential, interval.median_length)
        second_interval = Feature(interval.gene, interval.stack, interval.type, interval.order, 1, last(interval) - ctree.length, zero(Int8), interval.essential, interval.median_length)
        push!(ctree.ftree, first_interval)
        ctree.wrapped_intervals[first_interval] = interval
        push!(ctree.ftree, second_interval)
        ctree.wrapped_intervals[second_interval] = interval
    end
end

function Base.intersect(ftree::FeatureTree, btree::BlockTree)::Vector{Pair{Feature,AlignedBlock}}
    intervals = collect(intersect(ftree.ftree, btree.btree))
    return unique(map(x -> Pair(get(ftree.wrapped_intervals, x[1], x[1]), get(btree.wrapped_intervals, x[2], x[2])), intervals))
end

# entire set of Features for one strand of one genome
struct FeatureArray
    genome_id::String
    genome_length::Int32
    strand::Char
    feature_tree::FeatureTree
end

struct FeatureTemplate
    path::String  # similar to .sff path
    essential::Bool
    median_length::Float32 # median length of feature
end

# extended Feature struct to add prediction info
mutable struct SFF_Feature
    feature::Feature
    relative_length::Float32
    stackdepth::Float32
    gmatch::Float32
    feature_prob::Float32
    coding_prob::Float32
end

# alters reference_feature_count Dictionary
function read_sff_features!(file::String, reference_feature_counts::Dict{String,Int})::FwdRev{FeatureArray}
    open(file) do f
        header = split(readline(f), '\t')
        genome_id = header[1]
        genome_length = parse(Int32, header[2])
        r_features = FeatureTree(IntervalTree{Int32,Feature}(), genome_length, Dict{Feature,Feature}())
        f_features = FeatureTree(IntervalTree{Int32,Feature}(), genome_length, Dict{Feature,Feature}())
        while !eof(f)
            fields = split(readline(f), '\t')
            startswith(fields[1], "IR") && continue # don't use IR annotations
            startswith(fields[1], "unassigned") && continue # don't use unassigned annotations
            length(fields) ≥ 11 && length(fields[11]) > 0 && continue # don't use annotations with warnings
            feature = Feature(fields[1], parse(Int32, fields[3]), parse(Int32, fields[4]), parse(Int8, fields[5]))
            path = annotation_path(feature)
            count = get(reference_feature_counts, path, nothing)
            if isnothing(count)
                reference_feature_counts[path] = 1
            else
                reference_feature_counts[path] = count + 1
            end
            if fields[2][1] == '+'
                push!(f_features, feature)
            else
                push!(r_features, feature)
            end
        end
        f_strand_features = FeatureArray(genome_id, genome_length, '+', f_features)
        r_strand_features = FeatureArray(genome_id, genome_length, '-', r_features)
        return FwdRev(f_strand_features, r_strand_features)
    end
end

# part or all of a Feature annotated by alignment
struct Annotation
    genome_id::String
    path::String
    start::Int32
    length::Int32
    # offsets are the distance from the edge of the annotation to the end of 
    # the original feature
    offset5::Int32
    offset3::Int32
    # phase is the number of nucleotides to skip at the start of the sequence 
    # to be in the correct reading frame
    phase::Int8
end
datasize(a::Annotation) = sizeof(Annotation) + sizeof(a.path)# genome_id is shared + sizeof(a.genome_id)
datasize(v::Vector{Annotation}) = sum(datasize(a) for a in v)

function read_templates(file::String)::Dict{String,FeatureTemplate}
    if !isfile(file)
        error("\"$(file)\" is not a file")
    end
    if filesize(file) === 0
        error("no data in \"$(file)!\"")
    end
    templates = Dict{String,FeatureTemplate}()
    open(file) do f
        readline(f) # skip header
        for line in eachline(f)
            fields = split(line, '\t')
            template = FeatureTemplate(fields[1], parse(Bool, fields[2]), parse(Float32, fields[3]))
            templates[template.path] = template
        end
    end
    return templates
end

struct FeatureStack
    path::String # == template.path
    stack::CircularVector
    template::FeatureTemplate
end
datasize(f::FeatureStack) = sizeof(FeatureStack) + sizeof(f.path) + length(f.stack) * sizeof(Int8)

DFeatureStack = Dict{String,FeatureStack}
AFeatureStack = Vector{FeatureStack}

mutable struct SFF_Model
    gene::String
    gene_prob::Float32 # mean of probs of component features
    strand::Char
    gene_count::Int8
    exon_count::Int8
    features::Vector{SFF_Feature}
    warnings::Vector{String}
end

gene_length(model::Vector{Feature}) = begin
    maximum([f.start + f.length for f in model]) - minimum([f.start for f in model])
end

gene_length(model::SFF_Model) = begin
    maximum([m.feature.start + m.feature.length for m in model]) - minimum([m.feature.start for m in model])
end

function reverse_complement(r::UnitRange, glength::Int32)::UnitRange
    return range(mod1(glength - r.stop + 1, glength), length=length(r))
end

MaybeIR = Union{AlignedBlock,Nothing}


