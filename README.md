<img align="right" alt="Chloe" src="assets/logo-chloe-black.png">

# Chloë: Organelle Annotator

Chloë is optimised for annotating flowering plant (angiosperm) chloroplast genomes. If you would like to use Chloë to annotate chloroplast genomes from other plants (e.g. gymnosperms, ferns, lycophytes or bryophytes), please contact Ian Small (ian.small@uwa.edu.au) for access to future versions of Chloë.

This annotator is available online at: [https://chloe.plastid.org](https://chloe.plastid.org)

--- 

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Julia Projects](#julia-projects) 
- [Output Formats](#output-formats)
- [Authors](#authors)

## Installation

For installing Julia please follow link to [`juliaup`](https://julialang.org/downloads/).

To install the Chloë code as a local folder on your computer: 
```bash
git clone https://github.com/ian-small/chloe
```

#### Chloë References
Chloë references are required and can be cloned from the git repository into the same location as the `chloe` folder. If you want to save them somewhere else you need to specify the path to them.

```bash
git clone https://github.com/ian-small/chloe_references
```

The `Project.toml` file lists all the project
dependencies. To instantiate go into the `chloe` directory and type: 
```bash
julia --project=.
```
In the Julia REPL type `]` and `instantiate` at the julia prompt to install all the required
packages. 

## Usage

You can run Chloë from the terminal. To access the annotator help manual use:

```bash
julia --project=. chloe.jl annotate --help
```

Equivalently you can invoke Chloe with:

```bash
julia --project=. -e 'using Chloe; chloe_main()' annotate --help
```


For annotating single sequences (e.g. the test genome `NC_020019.1.fa` available in the folder `testfa` with the default output in `.sff` format:
```bash
julia --project=. chloe.jl annotate testfa/NC_020019.1.fa
```

For annotating all fasta file in a directory ending with `.fa` specifying the `.sff` output format: 

```bash
julia --project=. chloe.jl annotate --sff testfa/*.fa
```

This will create `.chloe.sff` files for each fasta file and write them back into the directory where the annotated fasta files are located.

To see what other commands are available:

```bash
julia --project=. chloe.jl --help
```

## Julia Projects
You can install Chloë as a Julia package and environment from within the Julia REPL. To create a project in your directory `myproject` initiate a Julia project and add Chloë as a package:

```bash
# create a Julia project in directory myproject
julia -e 'using Pkg; Pkg.generate("myproject")'
cd myproject
# add Chloe to the project
julia --project=. -e 'using Pkg; Pkg.add(url="https://github.com/ian-small/Chloe.jl.git")'
```

### Using the Julia REPL

Now you can start the Julia REPL and import the Chloë package.
As an example of how to annotate a single FASTA file that is in your project directory:
```julia
#start the julia REPL from the terminal in your projects folder then import Chloe
import Chloe
references = Chloe.ReferenceDb("cp")
outfile, uid = Chloe.annotate(references,  "NC_011032.1.fa") #run annotation on a file called NC_011032.1.fa located in your project folder
println(outfile) #print output in REPL
```

Or if you prefer you can use the commandline interface from the REPL to invoke Chloe:

```julia
import Chloe: chloe_main
chloe_main(["annotate", "-r", "nr", "nrdna.fa"])
```

--- 
For more recipes using Chloë see our [Recipes](https://github.com/ian-small/chloe/blob/master/RECIPES.md).


## Output formats

Internally, Chloë numbers each strand independently from its 5' end, and tracks features by (start, length) rather then by (start, stop). This avoids most of the issues with features crossing the arbitrary end of a circular genome. The `--sff` output of Chloë (`.sff` files) uses these conventions. For example, here's the start of a typical `.sff` output file:


<img src="assets/sff.png" width="600">


The header line gives the sequence name, the length in nucleotides, and the mean alignment coverage with the reference genomes.
Subsequent lines give information on a single feature or sub-feature.
The first column is a unique identifier, composed as follows:
gene name/gene copy (so if 2 or higher is a duplicate of another gene)/feature type/feature order (can be used to sort exons and introns into the correct order, even for transpliced genes)
Subsequent columns are: strand, start, length, phase;
Then 5 columns of interest if you want to understand why Chloë has predicted this particular feature: length relative to feature template, proportion of references that match, mean coverage of aligned genomes (out of 100), feature probability (from XGBoost model), coding probability (from XGBoost model)

The default output is GFF:

<img src="assets/gff.png" width="700">


By default, Chloë filters out features which are detected to have one of a set of problematic issues, or which have a feature probability of < 0.5.
You can retain these putative features by lowering the sensitivity threshold and asking for no filtering. For example, `chloe.jl annotate --sff --sensitivity 0 --no-filter` will retain all the features that Chloë was able to detect, including those that fail the checks. Features with issues will be flagged as warnings during the annotation:
```[ Warning: rps16/1 lacks a start codon
[ Warning: rps16/1 has a premature stop codon
[ Warning: rps16/1 CDS is not divisible by 3
```
and in the `.sff` output. Currently `--no-filter` has no effect if the `--sff` flag is not set.
<!-- Additional text to prevent interpretation as a header -->
---

### Authors

* Ian Small: ian.small@uwa.edu.au
* Ian Castleden: ian.castleden@uwa.edu.au
* Conny Hooper: cornelia.hooper@uwa.edu.au
