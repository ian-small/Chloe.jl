
# Development Notes

## Installing dependencies

The `Project.toml` file lists all the project
dependencies. From within the `chloe` directory, type `julia --project=.`
Then type `]instantiate` at the julia prompt to install all the required
packages.

## Chloë: Organelle Annotator

To run the annotator type:

```bash
julia --project=. chloe.jl annotate --help
```

For example:

```bash
julia --project=. chloe.jl annotate --gff testfa/*.fa
```

will create `.gff` files in the *current* directory.

This annotator is available online at: [https://chloe.plastid.org](https://chloe.plastid.org)

To see what other commands are available:

```bash
julia --project=. chloe.jl --help
```
## Chloë: Output formats

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

Most users will probably want to use `chloe.jl annotate -gff` to obtain the output in standard `.gff` format. 

By default, Chloe filters out features which are detected to have one of a set of problematic issues, or which have a feature probability of < 0.5.
You can retain these putative features by lowering the sensitivity threshold and asking for no filtering. For example, `chloe.jl annotate -s 0 --nofilter` will retain all the features that Chloe was able to detect, including those that fail the checks. Features with issues will be flagged as warnings during the annotation:
```[ Warning: rps16/1 lacks a start codon
[ Warning: rps16/1 has a premature stop codon
[ Warning: rps16/1 CDS is not divisible by 3
```
and in the `.sff` output. Currently `--nofilter` has no effect if the `--gff` flag is also set.

## Multithreading

Chloe will take advantage of multiple threads if possible. To benefit from this substantial speedup, specify the number of threads to use when starting Julia.
Using multiple threads is generally much faster than using multiple distributed processes (see the 'Distributed' section below).

For example:

```bash
julia --threads 4 --project=. chloe.jl annotate --gff testfa/*.fa
```
or

```bash
julia --threads auto --project=. chloe.jl annotate --gff testfa/*.fa
```

## Chloe as a Julia package

You can install Chloe as a Julia package.
Start Julia and type:

```julia
using Pkg;
Pkg.develop("{path/to/chloe/repo/directory}")
# *OR* directly from github
Pkg.add(url="https://github.com/ian-small/Chloe.jl.git")
```

This will make an entry for Chloë in the Manifest for Julia.
Now get Julia to compile it by typing `import Chloe` at the *julia* prompt.

You can easily remove Chloë as a package with:

```julia
using Pkg;
Pkg.rm("Chloe")
```

Installing Chloë as a (local) package allows you to take
advantage of Julia's precompilation.

## Distributed

* [Distributed](https://docs.julialang.org/en/v1/stdlib/Distributed/index.html)

```julia
using Distributed
import Chloe
addprocs(4)
@everywhere begin 
    using Chloe
    db = ReferenceDb("cp")
end
r = fetch(@spawnat :any annotate(db, "testfa/NC_020019.1.fa"))
# etc...
```

---

### Authors

* Ian Small: ian.small@uwa.edu.au
* Ian Castleden: ian.castleden@uwa.edu.au
