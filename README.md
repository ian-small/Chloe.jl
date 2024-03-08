
<img align="right" alt="Chloe" src="assets/logo-chloe-black.png">

# Chloë: Organelle Annotator

## Description

Chloë is optimised for annotating flowering plant (angiosperm) genomes. If you would like to use Chloë to annotate chloroplast genomes from other plants (e.g. gymnosperms, ferns, lycophytes or bryophytes), please contact Ian Small (ian.small@uwa.edu.au) for access to future versions of Chloë.

To use the web service please visit This annotator is available online at: [https://chloe.plastid.org](https://chloe.plastid.org)

--- 

## Table of Contents
- [Installing dependencies](#installing-dependencies)
    - [Git Clone Installation](#git-clone-installation)
    - [As Julia Package Installation](#as-julia-package-installation)  
- [Usage](#usage)
- [Output Formats](#output-formats)
- [Developer recipes](#developer-recipes)
    - [Multithreading](#multithreading)
    - [Distributed](#distributed)
    - [Chloë server](#chloe-server)
    - [Running Remotely](#running-remotely)
- [Authors] (#authors)


## Installing dependencies

### Git Clone installation

To install the Chloë code as a local folder:
```bash
git clone https://github.com/ian-small/chloe_references
```

#### Chloë References
Chloë references can be cloned from the repository to a folder in the same location as the Chloë folder if downloaded:
```bash
git clone https://github.com/ian-small/chloe_references
```

The `Project.toml` file lists all the project
dependencies. From within the `chloe` directory, type `julia --project=.`
Then type `]instantiate` at the julia prompt to install all the required
packages.


### As A Julia Package Installation
You can install Chloe as a Julia package from within the Julia REPL as either a environment package or a local package. If you just want to use the Chloe package as it is from its GitHub repository, you would use `add https://github.com/ian-small/chloe.git`. If you want to actively develop or modify the Chloe package, you might use `]dev path/to/chloe/` and work directly with the local source code (forked from the git repo).

How to: 
Start Julia and type `]` to get the package manager prompt.

Julia managed package:
```julia
add https://github.com/ian-small/chloe.git
```
Local package:
If you install it as a local package you need to git clone the repository to a location on your computer and then point Julia to this path using:
```julia
]dev {path/to/chloe/repo/directory}
```
This will make an entry for Chloë in the Manifest for Julia. Now get Julia to compile it by typing `import Chloe` at the *julia* prompt. 

---
You can easily remove Chloë as a package with:

```julia
]rm Chloe
```

Installing Chloë as a (local) package allows you to take
advantage of Julia's precompilation.


## Usage

To run the annotator type:

```bash
julia --project=. chloe.jl annotate --help
```

For example:

```bash
julia --project=. chloe.jl annotate -g testfa/*.fa
```

will create `.gff` files in the current directory.

To see what other commands are available:

```bash
julia --project=. chloe.jl --help
```
## Output formats

Internally, Chloe numbers each strand independently from its 5' end, and tracks features by (start, length)
rather then by (start, stop). This avoids most of the issues with features crossing the arbitrary end of a circular genome.
The default output of Chloe (`.sff` files) uses these conventions. For example, here's the start of a typical `.sff` output file:
`NC_020431.1	151328	78.914`
`accD/1/CDS/1	+	56703	1485	0	1.01	0.743	79.7	0.999	0.996`
The header line gives the sequence name, the length in nucleotides, and the mean alignment coverage with the reference genomes.
Subsequent lines give information on a single feature or sub-feature.
The first column is a unique identifier, composed as follows:
gene name/gene copy (so if 2 or higher is a duplicate of another gene)/feature type/feature order (can be used to sort exons and introns into the correct order, even for transpliced genes)
Subsequent columns are: strand, start, length, phase;
Then 5 columns of interest if you want to understand why Chloe has predicted this particular feature: length relative to feature template, proportion of references that match, mean coverage of aligned genomes (out of 100), feature probability (from XGBoost model), coding probability (from XGBoost model)

Most users will probably want to use `chloe.jl annotate -g` to obtain the output in standard `.gff` format. 

By default, Chloe filters out features which are detected to have one of a set of problematic issues, or which have a feature probability of < 0.5.
You can retain these putative features by lowering the sensitivity threshold and asking for no filtering. For example, `chloe.jl annotate -s 0 --nofilter` will retain all the features that Chloe was able to detect, including those that fail the checks. Features with issues will be flagged as warnings during the annotation:
```[ Warning: rps16/1 lacks a start codon
[ Warning: rps16/1 has a premature stop codon
[ Warning: rps16/1 CDS is not divisible by 3
```
and in the `.sff` output. Currently `--nofilter` has no effect if the `-g` flag is also set.


## Developer Recipes
### Multithreading

Chloe will take advantage of multiple threads if possible. To benefit from this substantial speedup, specify the number of threads to use when starting Julia.
Using multiple threads is generally much faster than using multiple distributed processes (see the 'Distributed' section below).

For example:

```bash
julia --threads 4 --project=. chloe.jl annotate -g testfa/*.fa
```
or

```bash
julia --threads auto --project=. chloe.jl annotate -g testfa/*.fa
```

### Distributed

Chloe has *[Distributed](https://docs.julialang.org/en/v1/stdlib/Distributed/index.html) capabilities meaning it can launch worker processes via the specified cluster manager. To use use julia's Distributed package, start julia (1.6) with 3 workers and load code:

```bash
julia --project=. -t 8  -p 3 -L src/dist/remote.jl
```
Now you can type:

```julia
using Distributed
# just read reference Data on remote workers

# note that REFS is not defined locally in the REPL!
reference_directory = "/some/directory"
r = @spawnat :any annotate_one(reference_directory, "testfa/NC_020019.1.fa")
io, uid = fetch(r)
sff = String(take!(io))
# this works too.., just tell Chloe the filename
r = @spawnat :any annotate_one(reference_directory, "testfa/NC_020019.1.fa", "write_to_this_file.sff")
```

This also works:

```julia
using Distrbuted
addprocs(3)

fasta = IOBuffer(read("testfa/NC_020019.1.fa", String))
io, uid = fetch(@spawnat :any annotate_one(reference_directory, fasta))
# get chloe sff as a string
sff = String(take!(io))
# *OR*
sff_filename, uid = fetch(@spawnat :any annotate_one(reference_directory, "testfa/NC_020019.1.fa", nothing))
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
end
# Note that neither REFS nor annotate_one is defined in the REPL
# ...but all is still good.
r = @spawnat :any annotate_one(reference_directory, "testfa/NC_020019.1.fa")
# etc...
```

This takes advantage of the precompilation of julia packages.
Also you don't need to be in the repo directory!

### Chloë Server

Running the chloe server. In a terminal type:

```bash
julia -t 8 --project=. distributed.jl --level=info --workers=4 \
     --broker="default"
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

### Running Remotely

The Chloë server can be run remotely through a ssh tunnel.

On the remote server:
`git clone ...` the chloe github repo and download the julia runtime (natch!).
*And* install all chloe package dependencies *globally* (see above).

Then -- on your puny laptop -- you can run something like:

```sh
ssh  you@bigserver -t -o ExitOnForwardFailure=yes -L 9476:127.0.0.1:9467 \
    'cd /path/to/chloe;
    JULIA_NUM_THREADS={BIGNUM} /path/to/bin/julia --project=. --startup-file=no --color=yes
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

### Authors

* Ian Small: ian.small@uwa.edu.au
* Ian Castleden: ian.castleden@uwa.edu.au
* Conny Hooper: cornelia.hooper@uwa.edu.au
