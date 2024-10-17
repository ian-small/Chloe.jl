include("globals.jl")

const KWARGS = ["no_transform", "sensitivity", "asgff3", "no_filter", "reference"]

struct ChloeConfig
    no_transform::Bool
    sensitivity::Real
    asgff3::Bool
    no_filter::Bool
    reference::String # cp|nr
    function ChloeConfig(;
        no_transform=false,
        sensitivity=DEFAULT_SENSITIVITY,
        asgff3::Bool=false,
        no_filter::Bool=false,
        reference::String="cp"
    )
        return new(no_transform, sensitivity, asgff3, no_filter, reference)
    end

    # needs to be V <: Any since this is coming from a JSON blob
    function ChloeConfig(dict::Dict{String,V} where {V<:Any})
        return ChloeConfig(; Dict(Symbol(k) => v for (k, v) in dict if k in KWARGS)...)
    end
end

function Base.show(io::IO, c::ChloeConfig)
    print(
        io,
        "ChloeConfig[no_transform=$(c.no_transform), sensitivity=$(c.sensitivity), no_filter=$(c.no_filter), asgff3=$(c.asgff3)], ref=$(c.reference)]"
    )
end

struct FeatureTemplate
    path::String  # similar to .sff path
    essential::Bool
    median_length::Float32 # median length of feature
    reference_strand::Char # strand feature is expected to be on in standard configuration of the genome/contig
end

abstract type AbstractReferenceDb end

mutable struct ReferenceDb <: AbstractReferenceDb
    lock::ReentrantLock
    gsrefsdir::String
    template_file::String
    templates::Union{Nothing,Dict{String,FeatureTemplate}}
    # gsrefhashes::Union{Nothing,Dict{String,Vector{Int64}}}
end

function ReferenceDb(reference_dir="cp")::ReferenceDb
    gsrefsdir = reference_dir
    if reference_dir == "cp"
        gsrefsdir = normpath(joinpath(CHLOE_REFS_DIR, "cprefs"))
    elseif reference_dir == "nr"
        gsrefsdir = normpath(joinpath(CHLOE_REFS_DIR, "nrefs"))
    end
    template = normpath(joinpath(gsrefsdir, DEFAULT_TEMPLATE))
    verify_refs(gsrefsdir, template)
    return ReferenceDb(ReentrantLock(), gsrefsdir, template, nothing)
end

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
            template = FeatureTemplate(fields[1], parse(Bool, fields[2]), parse(Float32, fields[3]), first(fields[4]))
            templates[template.path] = template
        end
    end
    return templates
end

function get_templates(db::ReferenceDb)
    lock(db.lock) do
        if isnothing(db.templates)
            db.templates = read_templates(db.template_file)
        end
        return db.templates
    end
end

function get_single_reference!(
    db::ReferenceDb,
    refID::AbstractString,
    reference_feature_counts::Dict{String,Int}
)::SingleReference
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

function verify_refs(gsrefsdir, template)
    # used by master process to check reference directory
    # *before* starting worker processes...
    if !isdir(gsrefsdir)
        msg = "Reference directory \"$(gsrefsdir)\" is not a directory!"
        @error msg
        throw(ArgumentError(msg))
    end
    if !isfile(template)
        msg = "template file \"$(template)\" does not exsit!"
        @error msg
        throw(ArgumentError(msg))
    end
end

# alters reference_feature_count Dictionary
function read_single_reference!(
    refdir::String,
    refID::AbstractString,
    reference_feature_counts::Dict{String,Int}
)::SingleReference
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
