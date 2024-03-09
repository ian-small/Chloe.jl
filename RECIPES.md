# Recipies

## Installing Julia

Follow the directions at the link to install [`juliaup`](https://julialang.org/downloads/).

## Creating a Julia Project

```sh
# create a Julia project in directory myproject
julia -e 'using Pkg; Pkg.generate("myproject")'
cd myproject
# add Chloe to the project
julia --project=. -e 'using Pkg; Pkg.add("https://github.com/ian-small/chloe.git")'
```
Get the Chloe database (This can be placed anywhere you want really)

```sh
git clone https://github.com/ian-small/chloe_references
# *OR* use julia
julia -e 'import Pkg; Pkg.GitTools.clone(stdout, "https://github.com/ian-small/chloe_references", "chloe_references")'
```

## CommandLine

Annotate fasta files from the command line.

```sh
# note the '--'
julia --project=. -e 'using Chloe; chloe_main()' -- \
    annotate --reference=/path/to/chloe_references *.fa
```
This will annotate files one-by-one.

## Using the Julia REPL

Annotate a single FASTA file.

```julia
import Chloe
# to quieten Chloe set the logging level:
# import Logging
# Logging.disable_logging(Logging.Info) # Disable debug and info

references = Chloe.ReferenceDbFromDir("/path/to/chloe_references")

outfile, uid = Chloe.annotate(references,  "NC_011032.1.fa")

println(outfile)
```

Write to buffer instead of to a file.

```julia
import Chloe
references = Chloe.ReferenceDbFromDir("/path/to/chloe_references")
io, uid = Chloe.annotate(references, "NC_011032.1.fa", nothing, IOBuffer())
# show .sff content
println(String(take!(io)))
```

Read from an already open fasta file.


```julia
import Chloe
references = Chloe.ReferenceDbFromDir("/path/to/chloe_references")
outfile, uid = open("NC_011032.1.fa", "r") do io
    Chloe.annotate(references, io)
end
```
## [Distributed](https://docs.julialang.org/en/v1/stdlib/Distributed/index.html)

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
    references = Chloe.ReferenceDbFromDir("/path/to/chloe_references")
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
    references = ReferenceDbFromDir("/path/to/chloe_references")
end
r = fetch(@spawnat :any annotate(references, "NC_011032.1.fa"))
println(r)
```

## Server

Run

```sh
julia -t8 --project=. -e 'using Chloe; distributed_main()' -- --level=info --workers=4 --broker=default --reference=/path/to/chloe_references
```

You can interact with this server using [JuliaWebAPI](https://github.com/JuliaWeb/JuliaWebAPI.jl)


```julia
using JuliaWebAPI
import Chloe
i = APIInvoker(Chloe.ZMQ_ENDPOINT)
apicall(i, "ping") # ping the server to see if is listening.
# annotate a file
ret = apicall(i, "annotate", read("NC_011032.1.fa",String))
code, data = ret["code"], ret["data"]
@assert code == 200
# actual filename written and total elapsed
# time in ms to annotate
sff_fname, elapsed_ms = data["filename"], data["elapsed"]

#stop the server....
apicall(i, "exit")
```
