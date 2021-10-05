

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

struct ChloeConfig
    numrefs::Int
    sensitivity::Real
    references::Union{Nothing,Vector{AbstractString}}
    to_gff3::Bool
    nofilter::Bool
end

function ReferenceDb(;refsdir="default",  hashfile="default",
    template="default")
    if refsdir == "default"
        refsdir = normpath(joinpath(REPO_DIR, "..", DEFAULT_REFS))
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
function get_minhashes(db::ReferenceDb, config::ChloeConfig)
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


const KWARGS = ["numrefs", "sensitivity", "to_gff3", "nofilter", "references"]

function ChloeConfig(;numrefs=DEFAULT_NUMREFS, sensitivity=DEFAULT_SENSITIVITY,
    references::Union{Nothing,Vector{AbstractString}, Vector{String}}=nothing, to_gff3::Bool=false, nofilter::Bool=false)
    return ChloeConfig(numrefs, sensitivity, references, to_gff3, nofilter)
end

# needs to be V <: Any since this is comming from a JSON blob
function ChloeConfig(dict::Dict{String,V} where V <: Any)
    function cvt(name, v)
        if name == "references"
            # if we actually have references then v
            # will be Any["str1","str2"]. convert to Vector{String}
            return string.(v)
        end
        # integers and float are ok
        v
    end
    return ChloeConfig(;Dict(Symbol(k) => cvt(k, v) for (k, v) in dict if k in KWARGS)...)
end
function Base.show(io::IO, c::ChloeConfig)
    print(io, "ChloeConfig[numrefs=$(c.numrefs), sensitivity=$(c.sensitivity), nofilter=$(c.nofilter), references=$(c.references)]")
end

function verify_refs(refsdir, hashfile, template)
    # used by master process to check reference directory
    # *before* starting worker processes...

    # TODO: read json file and really check...
    if !isdir(refsdir)
        msg = "Reference directory $(refsdir) is not a directory!"
        @error msg
        throw(ArgumentError(msg))
    end

    for f in [hashfile, template, joinpath(refsdir, "ReferenceOrganisms.json")]
        if !isfile(f)
            msg = "missing file: $f"
            @error msg
            throw(ArgumentError(msg))
        end
    end
    files = readdir(refsdir)    
    sff = findall(x -> endswith(x, ".sff"), files)
    if length(sff) == 0
        msg = "no reference .sff files in $(refsdir)!"
        @error msg
        throw(ArgumentError(msg))
    end
end

function read_single_reference!(refdir::String, refID::AbstractString, reference_feature_counts::Dict{String,Int})::SingleReference
    if !isdir(refdir); refdir = dirname(refdir); end
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
        ref_features = read_features!(normpath(joinpath(refdir, refID * ".sff")), reference_feature_counts)
        SingleReference(refID, CircularSequence(FASTA.sequence(ref)), ref_features)
    end
end