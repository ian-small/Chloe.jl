# Chloe: Organelle Annotator

To run the annotator type:

```
julia annotate.jl --help
```

## Chloe Server

Running the chloe server. In a terminal type:

```bash
JULIA_NUM_THREADS=4 julia src/chloe_svr.jl --level=info
```

In another terminal start julia:

```julia
using JuliaWebAPI

i = APIInvoker("tcp://127.0.0.1:9999")
# fasta and output should be relative to the server's working directory!
ret = apicall(i, "chloe", fastafile, outputfile) # outputfile is optional
code, data = ret["code"], ret["data"]
@assert code === 200
fname, elapsed_ms = data["filename"], data["elapsed"]
# to terminate the server
apicall(i, ":terminate")
```

## Installing depenencies

Start julia -- in this directory -- and type `]` then type:

```
pkg> activate .
pkg> instantiate
pkg> status
```


### Notes:

See:

* http://zguide.zeromq.org/py:all#Multithreading-with-ZeroMQ

Possible useful REPL packages

* add Revise
* add OhMyREPL
* `@code_warntype f()` check type system

