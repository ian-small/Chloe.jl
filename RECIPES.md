# Recipies

## Creating a Julia Project

```bash
# create a Julia project in directory myproject
julia -e 'using Pkg; Pkg.generate("myproject")'
cd myproject
# add Chloe to the project
julia --project=. -e 'using Pkg; Pkg.add("https://github.com/ian-small/chloe.git")'
# get Chloe database (This can be placed anywhere you want really)
git clone https://github.com/ian-small/chloe_references
```

## CommandLine

Annotate a fasta file from the command line.

```bash
# note the '--'
julia --project=. -e 'using Chloe; chloe_main()' -- annotate NC_011032.1.fa
```

## Simple

Annotate a single FASTA file.

```julia
import Chloe
# to quieten Chloe set the logging level:
# import Logging
# Logging.disable_logging(Logging.Info) # Disable debug and info

references = Chloe.ReferenceDbFromDir("chloe_references")

outfile, uid = Chloe.annotate(references,  "NC_011032.1.fa")

println(outfile)
```

Write to buffer instead of to a file.

```julia
import Chloe
references = Chloe.ReferenceDbFromDir("chloe_references")
io, uid = Chloe.annotate(references, "NC_011032.1.fa", nothing, IOBuffer())
# show .sff content
println(String(take!(io)))
```

Read from an already open fasta file.


```julia
import Chloe
references = Chloe.ReferenceDbFromDir("chloe_references")
outfile, uid = open("NC_011032.1.fa", "r") do io
    Chloe.annotate(references, io)
end
```
## Distributed

It's easy to annotate multiple fasta files in parallel

```julia
using Distributed
# add workers
addprocs(2)

@everywhere begin
    import Chloe
    # to quieten Chloe set the logging level:
    # import Logging
    # Logging.disable_logging(Logging.Info) # Disable debug and info
    references = Chloe.ReferenceDbFromDir("chloe_references")
end

fasta_directory = "fastas"

pmap = f -> l -> map(f, l)
# find all fasta files in a directory
fastas = readdir(dir) |> filter(f -> endswith(f, r"\.fa")) |> pmap(f -> joinpath(fasta_directory, f))

# outputs is the list of output files
outputs = @distributed (vcat) for fasta = fastas
    # note that `references` in on the worker process
    output, uid = Chloe.annotate(references, fasta, nothing, fasta * ".sff")
    [output]
end
```

```julia
using Distributed
addprocs(4)
@everywhere  begin
    using Chloe
    references = ReferenceDbFromDir("~/chloe_references")
end
r = fetch(@spawnat :any annotate(references, "NC_011032.1.fa"))
println(r)
```
