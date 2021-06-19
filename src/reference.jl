

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
    lock(db.lock) do
        if isnothing(db.templates)
            db.templates = read_templates(db.template_file)
        end
        return db.templates
    end
end
function get_minhashes(db::ReferenceDb)
    lock(db.lock) do
        if isnothing(db.templates)
            db.refhashes = readminhashes(db.hash_file)
        end
        return db.refhashes
    end
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

const KWARGS = ["numrefs", "sensitivity", "to_gff3", "nofilter"]

function ChloeConfig(;numrefs=DEFAULT_NUMREFS, sensitivity=DEFAULT_SENSITIVITY, to_gff3::Bool=false, nofilter::Bool=false)
    return ChloeConfig(numrefs, sensitivity, to_gff3, nofilter)
end

# needs to be V <: Any since this is comming from a JSON blob
function ChloeConfig(dict::Dict{String,V} where V <: Any)
    return ChloeConfig(;Dict(Symbol(k) => v for (k, v) in dict if k in KWARGS)...)
end
function Base.show(io::IO, c::ChloeConfig)
    print(io, "ChloeConfig[numrefs=$(c.numrefs), sensitivity=$(c.sensitivity), nofilter=$(c.nofilter)]")
end