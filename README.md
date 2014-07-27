# How to use

This code is a bare-bones implementation of the Gillespie SSA algorithm (on a uniform time grid). The simulation code lives in the module directory, and is run from the IJulia notebook. To use it, you need to:

- install [julia](http://julialang.org/downloads/) (version 0.3.0 or greater).  
- install ipython notebook (installing [anaconda](http://docs.continuum.io/anaconda/install.html) is the easiest way).
- add the ijulia package (`Pkg.add("IJulia")` in the julia prompt)

The example system is set up to stochastically simulate a cell that grows exponentially in legth as it exponentially 
increases the amount of a ceratin kind of protein it contains (e.g.  gene expression in e. coli where the gene 
facilitates its own production).  The cell undergoes (noisily) periodic cell division based on its length, and at cell division, the amount of protein the cell 
retains is stocasticallly determined.  The parameters that govern these rates and noise can be changed in the ijulia 
notebook.

The system being simulated can also be changed by altering the initial state vector and the stoichiometry matrix to 
appropriately model any desired stochastic chemical kintetics.  For example, instead of the reduced model for 
constituative growth above, we could model the central dogma more explicitly by including mRNA in the process.
