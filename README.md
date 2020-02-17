# A Julia interface to PGPlot


[`PGPlot`](http://www.astro.caltech.edu/~tjp/pgplot/) is a library written by
Tim Pearson for scientific plots.  The `PGPlot.jl` package provides a Julia
interface to this library.
If you are tired by the slowness of other Julia plotting packages, you may give
a try to this one!

There are two layers a low-level interface and a higher level one.  The low
level interface is accessible via:

```using
using PGPlot.Bindings
```

and provides bindings for all PGPlot routines, although to a somewhat higher
level than the original FORTRAN or C routines.  You can see examples
[here](test/demo.jl).


The high level interface is a work in progress.


## Installation

For now installation is by-hand (it will be automated in a near future).

1. Make sure PGPlot (package `pgplot5` on Ubuntu and derived systems) is
   properly installed.  There is a Giza library which aims at providing a
   drop-in replacement of PGPlot using Cairo for anti-aliased graphics but all
   tests have been done with the original PGPlot.

2. Clone the `PGPlot.jl` repository.  From Julia interpreter, hit the ] key to
   switch to the package manager REPL (you should get a `... pkg>` prompt) and
   type:

   ```julia
   ... pkg> add https://github.com/emmt/PGPlot.jl
   ```

3. Edit the `deps.jl` file in the package `deps` sub-directory to reflect the
   location of the dynamic PGPlot library.

4. Try the numerous demos [here](test/demo.jl).
