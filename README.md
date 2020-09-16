
<img align="right" alt="Chloe" src="assets/logo-chloe-black.png">


# Chloë: Organelle Annotator

To run the annotator or write gff3 or create suffix array files type:

```bash
julia chloe.jl --help
# or for a specific command e.g.
julia chloe.jl annotate --help
```
Of if you have added Chloe as a package: `julia -e 'using Chloe; cmd_main()' -- annotate --help`

(See installing dependencies below)

For example:

```bash
julia chloe.jl annotate testfa/*.fa
```

Will create `.sff` files in the testfa directory.

This annotator is available online at: https://chloe.plantenergy.edu.au

## Installing dependencies

There is a `Project.toml` file that contains all the project
dependencies.

Start julia and type `]` to get the package manager prompt. Then type:

```
(@v1.5) pkg> dev {path/to/chloe/repo/directory}
```
This will make an entry for Chloe in the Manifest for julia. More
importantly it will install all the dependencies.
Now get julia to compile it by typing `import Chloe` at the *julia* prompt.

You can remove Chloe with:

```
(@v1.5) pkg> rm Chloe
```


## Chloë Server

Running the chloe server. In a terminal type:

```bash
JULIA_NUM_THREADS=8 julia distributed.jl --level=info --workers=4 \
     --broker=ipc:///tmp/chloe-client
```
(Julia as of 1.4 refuses to use more threads that the number of CPUs on your machine:
`Sys.CPU_THREADS` or `python -c 'import multiprocessing as m; print(m.cpu_count())'`)

In another terminal start julia:

```julia
using JuliaWebAPI

i = APIInvoker("ipc:///tmp/chloe-client");
apicall(i, "ping") # ping the server to see if is listening.

# fasta and output should be relative to the server'
# working directory, or specify absolute path names! yes "chloe" 
# should be "annotate" but...
ret = apicall(i, "chloe", fastafile, outputfile) # outputfile is optional
code, data = ret["code"], ret["data"]
@assert code === 200
# actual filename written and total elapsed
# time in ms to annotate
sff_fname, elapsed_ms = data["filename"], data["elapsed"]
# to terminate the server cleanly (after finishing any work)
apicall(i, "exit")

```

The *actual* production configuration uses `distributed.jl` 
(for threading issues) and runs
the server as a client of a DEALER/ROUTER server
(see `bin/broker.py` or `src/broker.jl` and the `Makefile`). It *connects* to the
DEALER end on `tcp://127.0.0.1:9467`. The
[chloe website](https://chloe.plantenergy.edu.au)
connects to `ipc:///tmp/chloe-client` which
is the ROUTER end of broker. In this setup
you can run multiple chloe servers connecting
to the same DEALER.

**Update**: you can now run a broker with julia as `julia src/broker.jl`
*or* specify `--broker=URL` to `distrbuted.jl`. No
python required.

The worker process can be made to share the reference Data using memory mapped data files.
You can create these by running:

```sh
julia chloe.jl mmap reference_1116/*.fa
```


## Distributed

* https://docs.julialang.org/en/v1/stdlib/Distributed/index.html

You can of course use julia's Distributed package.

Start julia with 3 workers and load code:

`JULIA_NUM_THREADS=8 julia -p 3 -L src/remote.jl`

Now you can type:

```julia
using Distributed
# just read reference Data on remote workers
@everywhere workers() begin
    global REFS = readReferences("reference_1116", "optimised_templates.v2.tsv")
end
# get a fasta file
fasta = IOBuffer(read("testfa/NC_020019.1.fa", String))
# note that REFS is not defined locally!
r = @spawnat :any annotate_one(REFS, fasta)
io, uid = fetch(r)
sff = String(take!(io))
# this works too.., just tell Chloe the filename 
r = @spawnat :any annotate_one(REFS, "testfa/NC_020019.1.fa")
r = @spawnat :any annotate_one(REFS, "testfa/NC_020019.1.fa", "write_to_this_file.sff")
```

This also works:

```julia
using Distrbuted
addprocs(3)
@everywhere workers() begin
    include("src/remote.jl")
    REFS = readReferences("reference_1116", "optimised_templates.v2.tsv")
end
fasta = IOBuffer(read("testfa/NC_020019.1.fa", String))
io, uid = fetch(@spawnat :any annotate_one(REFS, fasta))
# get chloe sff as a string
sff = String(take!(io))
# *OR*
sff_filename, uid = fetch(@spawnat :any annotate_one(REFS, "testfa/NC_020019.1.fa", nothing))
# sff_filename is where chloe wrote the data:
# in this case NC_020019.1.sff in the local directory
# instead of `nothing` specify an actual filename.
```

The

## Running Remotely

The Chloë server can be run remotely through a ssh tunnel.

On the remote server:
`git clone ...` the chloe github repo and download the julia runtime (natch!).
*And* install all chloe package dependencies *globally* (see above).

Then -- on your puny laptop -- you can run something like:

```sh
ssh  you@bigserver -t -o ExitOnForwardFailure=yes -L 9476:127.0.0.1:9467 \
    'cd /path/to/chloe; JULIA_NUM_THREADS={BIGNUM} /path/to/bin/julia --startup-file=no --color=yes distributed.jl  
    --broker=tcp://127.0.0.1:9467 -l info --workers=4'
```
The port `9467` is an entirely random (but hopefully unused both on
the remote server and locally) port number. The broker port *must* match
the ssh port specified by `-L`. `{BIGNUM}` is the enormous number
of CPUs your server has ;).

Since the remote server has no access to the local filesystem you need
to use `annotate` instead of `chloe` to annotate your your
fasta files e.g:

```julia
using JuliaWebAPI
i = APIInvoker("tcp://127.0.0.1:9467")
# read in the entire fasta file
fasta = read("testfa/NC_020019.1.fa", String)
ret = apicall(i, "annotate", fasta)
code, data = ret["code"], ret["data"]
@assert code === 200
sff = data["sff"] # sff file as a string
# terminate the server
apicall(i, "exit")
```

### Developer Notes:

To stop julia vomiting unhelpful stacktraces when `^Ctrl-C`ing 
run julia with `--handle-signals=no`. Don't know what it does
but `distributed.jl` will just exit on Ctrl-C.

But don't send a `kill -INT` this will not clean up the background
broker (if it's running)

See:

* http://zguide.zeromq.org/py:all#Multithreading-with-ZeroMQ

Possibly useful REPL packages

* add Revise: reload edited files within REPL
* add OhMyREPL: pretty print code
* `@code_warntype f()` check type system
* add ProfileView: https://github.com/timholy/ProfileView.jl


### Authors

* Ian Small: ian.small@uwa.edu.au
* Ian Castleden: ian.castleden@uwa.edu.au
