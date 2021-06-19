

# export ReferenceDb, get_minhashes, get_templates,get_single_reference!
include("globals.jl")
mutable struct ReferenceDb
    lock::ReentrantLock
    refsdir::String
    template_file::String
    hash_file::String
    templates::Union{Nothing,Dict{String,FeatureTemplate}}
    refhashes::Union{Nothing,Dict{String,Vector{Int64}}}
end

function ReferenceDb(;refsdir="default",  hashfile="default",
    template="default")
    if refsdir == "default"
        refsdir = normpath(joinpath(HERE, "..", DEFAULT_REFS))
    end
    if hashfile == "default"
        hashfile = normpath(joinpath(refsdir, DEFAULT_HASHES))
    end
    if template == "default"
        template = normpath(joinpath(refsdir, DEFAULT_TEMPLATE))
    end
    return ReferenceDb(ReentrantLock(), refsdir, template, hashfile, nothing, nothing)
end

function get_templates(db::ReferenceDb)
    if isnothing(db.templates)
        lock(db.lock) do
            db.templates = read_templates(db.template_file)
        end
    end
    db.templates
end
function get_minhashes(db::ReferenceDb)
    if isnothing(db.templates)
        lock(db.lock) do
            db.refhashes = readminhashes(db.hash_file)
        end
    end
    db.refhashes
end
function get_single_reference!(db::ReferenceDb, refID::AbstractString, reference_feature_counts::Dict{String,Int})::SingleReference
    read_single_reference!(db.refsdir, refID, reference_feature_counts)
end

struct ChloeConfig
    numrefs::Int
    sensitivity::Real
    to_gff3::Bool
    nofilter::Bool
end

function ChloeConfig(;numrefs=DEFAULT_NUMREFS, sensitivity=DEFAULT_SENSITIVITY, to_gff3::Bool=false, nofilter::Bool=false)
    return ChloeConfig(numrefs, sensitivity, to_gff3, nofilter)
end

function ChloeConfig(dict::Dict{String,Any})
    return ChloeConfig(;Dict(Symbol(k) => v for (k, v) in dict)...)
end