# CHANGES

## Sept 28 2020

* new `julia chloe.jl rotate` command to rotate genomes to standard position
* User API for alignments.

Subject to name changes you can get a bacic alignment with:

```julia
import Chloe
# you can use a mmapped file too!
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

* `compareSubStrings` is now ~50-400x faster! (most of the alignment
   time is spent in sorting all the (sub)sequences however).
* Reference data files can be memory mapped so all workers can now share reference data
  (use `julia chloe.jl mmap reference_1116/*.fa` to generate). **Big win here**.
  Means we can really *up* the number of reference genomes...
* Implemented an `AbstractString` called `MMappedString` to memory map
  a file as a julia String. (This really added to my gray hairs!) See below.
* Possibility of using only forward reference data (memory mapping makes this less of
  a win)
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

The other main benefit of this is that I've been forced to
create a nifty `AbstractString` implementation to memory map String Data
(without reading the data in up front). This was much harder than
I thought.

```julia

   import Chloe: MMappedString, ASCII
   import Mmap

   f = open("mysequencefile.txt")
   s = MMappedString(Mmap.map(f))
   # If you are *sure* that all the data in the file is ASCII use
   s = MMappedString{ASCII}(Mmap.map(f))
   # this is what gives me the 50x win noted above.
```

Now you have a nice String representation of the file
without reading a thing!

The annotation interface is basically

`annotate_one(REFERENCE, fasta)`

Stupidly I was reading the REFERENCE data in the master process
and copying it to worker processes on every call
(I thought julia Distributed was smarter... I was wrong). Now
each worker has its own copy of REFERENCE (see memory mmapping above)
so things are a little speedier.
