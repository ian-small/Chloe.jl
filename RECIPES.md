# Developer Recipies
- [Distributed](#distributed)
- [Server](#server)
- [Running Remotely](#running-remotely)




## [Distributed](https://docs.julialang.org/en/v1/stdlib/Distributed/index.html)

Chloë can use the Julia [Distributed](https://docs.julialang.org/en/v1/stdlib/Distributed/index.html) computing to annotate fasta files by distributing thw work to parallel processes. Steps are, you add worker processes, import Chloë with necessary references, then distribute the annotation process across workers, generating output files in parallel.

```julia
using Distributed
# add workers
addprocs(2)

@everywhere begin
    using Chloe
    # to quieten Chloe set the logging level:
    # import Logging
    # Logging.disable_logging(Logging.Info) # Disable debug and info
    references = ReferenceDb("cp")
end

fasta_directory = "fastas"

pmap = f -> l -> map(f, l)
# find all fasta files in a directory
fastas = readdir(fasta_directory) |> filter(f -> endswith(f, r"\.fa")) |> pmap(f -> joinpath(fasta_directory, f))

# outputs is the list of output files
outputs = @distributed (vcat) for fasta = fastas
    # note that `references` in on the worker process
    output, uid = annotate(references, fasta, nothing, fasta * ".sff")
    [output]
end
```

```julia
using Distributed
addprocs(4)
@everywhere  begin
    using Chloe
    references = ReferenceDb("cp")
end
r = fetch(@spawnat :any annotate(references, "NC_011032.1.fa"))
println(r)
```
--- 
### Authors

* Ian Small: ian.small@uwa.edu.au
* Ian Castleden: ian.castleden@uwa.edu.au
* Conny Hooper: cornelia.hooper@uwa.edu.au
