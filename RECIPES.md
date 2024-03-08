# Recipies

## Creating a Project

```bash
julia -e 'using Pkg; Pkg.generate("myproject")' \
    && cd myproject \
    && julia --project=. -e 'using Pkg; Pkg.add("https://github.com/ian-small/chloe.git")'
```

## Simple

Annotate a single FASTA file

```julia
import Chloe
import Logging

# Logging.disable_logging(Logging.Info) # Disable debug and info

# cd && git clone https://github.com/ian-small/chloe_references
db = Chloe.ReferenceDbFromDir("~/chloe_references")

outfile, uid = Chloe.annotate(db,  "NC_011032.1.fa")

println(outfile)
```
## Distributed

```julia
using Distributed
# add workers
addprocs(2)

@everywhere begin
    import Chloe
    # import Logging
    # Logging.disable_logging(Logging.Info) # Disable debug and info
    references = Chloe.ReferenceDbFromDir("~/chloe_references")
end

fasta_directory = "fastas"

pmap = f -> l -> map(f, l)
fastas = readdir(dir) |> filter(f -> endswith(f, r"\.fa")) |> pmap(f -> joinpath(fasta_directory, f))
outputs = @distributed (vcat) for fname = fastas
    # note that references in on a worker process
    output, uid = Chloe.annotate(references, fname, nothing, fname * ".sff")
    [output]
end
```

```julia
using Distributed
addprocs(4)
@everywhere  begin
    using Chloe
    db = ReferenceDbFromDir("~/chloe_references")
end
r = fetch(@spawnat :any annotate(db, "NC_011032.1.fa"))
println(r)
```
