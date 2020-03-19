# Chloe: Organelle Annotator

# Chloe Server

Running the chloe server

```bash
JULIA_NUM_THREADS=4 julia src/chloe_svr.jl --level=info
```

Invoking values on the server

```julia
using JuliaWebAPI

i = APIInvoker("tcp://127.0.0.1:9999")
# fasta and output should be relative to the server's working directory!
ret = apicall(i, "chloe", fasta, output)
@assert ret["code"] == 200
fname = ret["data"]["data"]
```