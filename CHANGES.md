# CHANGES

## Feb 4 2021

* new biojulia version with new alignment 'engine'
* should be faster and more accurate than before
* can pick the best references to use out of a large collection by comparing minhashes
* new `julia chloe.jl minhash [file or directory]` command to generate minhashes for a custom reference collection
* use `julia chloe.jl annotate -r [reference fasta files] --minhashes [reference minhashes]` to direct chloe to use your custom reference collection
* still to do: dynamic calculation of annotation stack thresholds
  for references other than the default set, the default templates are likely to be suboptimal

## Sept 28 2020

* new `julia chloe.jl rotate` command to rotate genomes to standard position
* User API for alignments.

Subject to name changes you can get a basic alignment with:

```julia
import Chloe
r1 = Chloe.createTargetReference("testfa/NC_020361.1.fa")
r2 = Chloe.createTargetReference("testfa/NC_020431.1.fa")
# Align the two sequences
# haven't exported this function to the top level...
blocks = Chloe.Annotator.align(r1,r2);
# average coverage
Chloe.Annotator.avg_coverage(r1, blocks) # uses blocks.forward
# or get counts directly
Chloe.Annotator.blockCoverage(blocks.reverse.forward)
# there is forward.forward,forward.reverse,reverse.forward and reverse.reverse
# inverted repeat
fwd, rev, len = Chloe.Annotator.inverted_repeat(r1)
```

## Sept 20 2020

### Executive Summary

* I've used an IntervalTree to allow quicker mapping of AlignedBlocks to
  reference Features (new dependency however).
* bugfix: `calc_maxlengths` was zipping the forward and reverse models viz:
  `zip(fstrand_models,rstrand_models)`. But there is no reason these two arrays
  are the same length. The change means I'm getting more `possible pseudogene` calls.
* *not* transfering 40MB of reference data for each annotation call. (see below).
* `lcps2AlignmentBlocks` was not checking the last elmement in lcps
* implemented an `iter_wrap` iterator to do proper wrapped searches for stop codons etc.

#### Gory Details

Chloe is a full julia package you can type
`]dev {path/to/chloe/repo}` to "install it". You get julia's
pre-compilation as a benefit.

The annotation interface is basically

`annotate_one(REFERENCE, fasta)`

Stupidly I was reading the REFERENCE data in the master process
and copying it to worker processes on every call
(I thought julia Distributed was smarter... I was wrong). Now
each worker has its own copy of REFERENCE (see memory mmapping above)
so things are a little speedier.
