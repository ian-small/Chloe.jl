include("globals.jl")
abstract type AbstractReferenceDb end

mutable struct ReferenceDb <: AbstractReferenceDb
    lock::ReentrantLock
    gsrefsdir::String
    template_file::String
    templates::Union{Nothing,Dict{String,FeatureTemplate}}
    gsrefhashes::Union{Nothing,Dict{String,Vector{Int64}}}
end

struct ChloeConfig
    numgsrefs::Int
    sensitivity::Real
    to_gff3::Bool
    nofilter::Bool
end

function ReferenceDb(; gsrefsdir="default", template="default")::ReferenceDb
    if gsrefsdir == "default"
        gsrefsdir = normpath(joinpath(REPO_DIR, "..", "..", DEFAULT_GSREFS))
    end
    if template == "default"
        template = normpath(joinpath(REPO_DIR, "..", "..", DEFAULT_TEMPLATE))
    end
    return ReferenceDb(ReentrantLock(), gsrefsdir, template, nothing, nothing)
end

function get_templates(db::ReferenceDb)
    lock(db.lock) do
        if isnothing(db.templates)
            db.templates = read_templates(db.template_file)
        end
        return db.templates
    end
end

function get_gsminhashes(db::ReferenceDb, config::ChloeConfig)
    config.numgsrefs < 1 && return nothing
    lock(db.lock) do
        if isnothing(db.gsrefhashes)
            db.gsrefhashes = readminhashes(normpath(joinpath(db.gsrefsdir, "reference_minhashes.hash")))
        end
        return db.gsrefhashes
    end
end


function get_single_reference!(db::ReferenceDb, refID::AbstractString, reference_feature_counts::Dict{String,Int})::SingleReference
    path = findfastafile(db.gsrefsdir, refID)
    if isnothing(path) || !isfile(path)
        msg = "unable to find $(refID) fasta file in $(db.gsrefsdir)!"
        @error msg
        throw(ArgumentError(msg))
    end
    open(path) do io
        ref = FASTA.Record()
        reader = FASTA.Reader(io)
        read!(reader, ref)
        sffpath = path[1:findlast('.', path)] * "sff" #assumes fasta files and sff files differ only by the file name extension
        if !isfile(sffpath)
            msg = "unable to find $(refID) sff file in $(db.gsrefsdir)!"
            @error msg
            throw(ArgumentError(msg))
        end
        ref_features = read_sff_features!(sffpath, reference_feature_counts)
        SingleReference(refID, CircularSequence(FASTA.sequence(LongDNA{4}, ref)), ref_features)

    end
end

const KWARGS = ["numgsrefs", "sensitivity", "to_gff3", "nofilter"]

function ChloeConfig(; numgsrefs=DEFAULT_NUMGSREFS, sensitivity=DEFAULT_SENSITIVITY,
    to_gff3::Bool=false, nofilter::Bool=false)
    return ChloeConfig(numgsrefs, sensitivity, to_gff3, nofilter)
end

# needs to be V <: Any since this is comming from a JSON blob
function ChloeConfig(dict::Dict{String,V} where {V<:Any})
    function cvt(name, v)
        if name == "references"
            # if we actually have references then v
            # will be Any["str1","str2"]. convert to Vector{String}
            return string.(v)
        end
        # integers and float are ok
        v
    end
    return ChloeConfig(; Dict(Symbol(k) => cvt(k, v) for (k, v) in dict if k in KWARGS)...)
end
function Base.show(io::IO, c::ChloeConfig)
    print(io, "ChloeConfig[numgsrefs=$(c.numgsrefs), sensitivity=$(c.sensitivity), nofilter=$(c.nofilter)]")
end

function verify_refs(gsrefsdir, template)
    # used by master process to check reference directory
    # *before* starting worker processes...
    if !isdir(gsrefsdir)
        msg = "Reference directory $(gsrefsdir) is not a directory!"
        @error msg
        throw(ArgumentError(msg))
    end
end

function read_single_reference!(refdir::String, refID::AbstractString, reference_feature_counts::Dict{String,Int})::SingleReference
    if !isdir(refdir)
        refdir = dirname(refdir)
    end
    path = findfastafile(refdir, refID)
    if isnothing(path)
        msg = "unable to find $(refID) fasta file in $(refdir)!"
        @error msg
        throw(ArgumentError(msg))
    end
    open(path) do io
        ref = FASTA.Record()
        reader = FASTA.Reader(io)
        read!(reader, ref)

        ref_features = read_sff_features!(normpath(joinpath(refdir, refID * ".sff")), reference_feature_counts)
        SingleReference(refID, CircularSequence(FASTA.sequence(LongDNA{4}, ref)), ref_features)
    end
end
