
<img align="right" alt="Chloe" src="assets/logo-chloe-black.png">

# Chloë: Organelle Annotator

To run the annotator or write gff3 or create suffix array files type:

```bash
julia chloe.jl --help
# or for a specific command e.g.
julia chloe.jl annotate --help
```

(See installing dependencies below)

For example:

```bash
julia chloe.jl annotate testfa/*.fa
```

Will create `.sff` files in the testfa directory.

This annotator is available online at: [https://chloe.plantenergy.edu.au](https://chloe.plantenergy.edu.au)

## Installing dependencies

There is a `Project.toml` file that contains all the project
dependencies.

To actually add these dependencies type
`julia bin/deps.jl`.

*or* run

```julia
import Pkg
(open("Project.toml") |> Pkg.TOML.parse)["deps"] |> keys |> collect |> Pkg.add
```

I really don't know why there isn't a command for this... :(

You can install Chloe as a julia
package too.
Start julia and type `]` to get the package manager prompt. Then type:

```julia
]dev {path/to/chloe/repo/directory}
```

This will make an entry for Chloë in the Manifest for julia.
Now get julia to compile it by typing `import Chloe` at the *julia* prompt.

You can easily remove Chloë as a package with:

```julia
]rm Chloe
```

Installing Chloë as a (local) package allows you to take
advantage of julia's precompilation.

## Distributed

* [Distributed](https://docs.julialang.org/en/v1/stdlib/Distributed/index.html)

You can of course use julia's Distributed package.

Start julia with 3 workers and load code:

`JULIA_NUM_THREADS=8 julia -p 3 -L src/remote.jl`

Now you can type:

```julia
using Distributed
# just read reference Data on remote workers
@everywhere workers() begin
    global REFS = readDefaultReferences()
end
# get a fasta file
fasta = IOBuffer(read("testfa/NC_020019.1.fa", String))
# note that REFS is not defined locally in the REPL!
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
    REFS = readDefaultReferences()
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

If you have installed Chloe as a (local) package the you can use:

```julia
using Distributed
addprocs(4)
@everywhere workers() begin
    using Chloe
    global REFS = readDefaultReferences()
end
# Note that neither REFS nor annotate_one is defined in the REPL
# ...but all is still good.
r = @spawnat :any annotate_one(REFS, "testfa/NC_020019.1.fa")
# etc...
```

This takes advantage of the precompilation of julia packages.
Also you don't need to be in the repo directory!

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
python required. (best to use `-b default` to select
this projects default endpoint (`ipc:///tmp/chloe-client`))

The worker process can be made to share the reference Data using memory mapped data files.
You can create these by running:

```sh
julia chloe.jl mmap reference_1116/*.fa
```

## Running Remotely

The Chloë server can be run remotely through a ssh tunnel.

On the remote server:
`git clone ...` the chloe github repo and download the julia runtime (natch!).
*And* install all chloe package dependencies *globally* (see above).

Then -- on your puny laptop -- you can run something like:

```sh
ssh  you@bigserver -t -o ExitOnForwardFailure=yes -L 9476:127.0.0.1:9467 \
    'cd /path/to/chloe;
    JULIA_NUM_THREADS={BIGNUM} /path/to/bin/julia --startup-file=no --color=yes
    distributed.jl --broker=tcp://127.0.0.1:9467 -l info --workers=4'
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

---

### Developer Notes

Nothing interesting beyond here....

Install package with

```julia
Pkg.clone("https://github.com:arabidopsis/chloe.git")
```

To stop julia vomiting unhelpful stacktraces when `^Ctrl-C`ing
run julia with `--handle-signals=no`. Don't know what it does
but `distributed.jl` will just exit on Ctrl-C.

But don't send a `kill -INT` this will not clean up the background
broker (if it's running)

See:

* [Multithreading in ZMQ](http://zguide.zeromq.org/py:all#Multithreading-with-ZeroMQ)

Possibly useful REPL packages

* add Revise: reload edited files within REPL
* add OhMyREPL: pretty print code
* `@code_warntype f()` check type system
* add ProfileView: [ProfileView.jl](https://github.com/timholy/ProfileView.jl)

from [stackoverflow](https://stackoverflow.com/questions/38825626/julia-transferring-methods-between-workers/39216340#39216340):

There is no way to send a subset of the methods in a package to another machine.
Very often methods refer to other types and functions in the same module, so the
system would have to at least send all dependencies as well. That could work, but
the bigger problem is deciding whose responsibility it is to distribute code,
and when. For example, initially your library might decide to send itself
(or parts of itself) to other nodes, but then the user might later want to do a
parallel map of your library functions, such that the whole library is needed on
every node. This gets very complex, so it is far simpler for everybody just to load
all needed code on all nodes as early as possible.

---

This only really is of interest with using the `add_worker` method that tries
to add new workers to the running server. If the server was
started by loading `Chloe` as a package then you can't add new workers by just sending
the required code: The new worker seems to be expecting a Chloe module.
Use `distributed.jl` if you want to expand workers dynamically.

### Authors

* Ian Small: ian.small@uwa.edu.au
* Ian Castleden: ian.castleden@uwa.edu.au
