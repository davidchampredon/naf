# naf
## New agent-based model for flu

Epidemiological model is coded in C++ and wrapped in R for output analysis. 

Preferred way to run the simulations:
* 0 [Optional]) Make sure the C++ code is correct by compiling it: `make naf`. Code can be run with `./naf`.
* 1) Build the R `naf` library: in folder `Rlibrary`, execute `./build_library`. If there is a problem at this stage only (if the previous optional stage is run without any error), the problem probably lies in `Rwrap.cpp`.
* 2 [Optional]) Fit the age distributions conditional on household size on the global age distribution: in folder `fit`, run `go-fit-hhszage`. Once the fit is done, generate the best fitted distributions by running the script `go-gen-ad-hhsz` in the folder `data/households`. The best fitted parameters must be copy-pasted in `gen-ad-hhsz.R` (TO DO: change this).
* 3) Fit the contact rate on arbitrary values of R0: edit and run the script `fit/go-fitR`. Several values for the baseline contact rate will be used to calculate R0. Pick the one that is close to the desired R0. (So, it's not an automatic fit).
* 4) Simulations with scenario comparison are launched on a local computer with the script `go-simul` in the `simul` folder.  
* 5) To launch multiple scenarios on large clusters (high-performance computing, HPC), execute in the following order:
 * read the check list `_checklist_HPC.txt`
 * Define the scenarios in `simul/scenario-builder.R`
 * Launch simulations with `simul/go-hpc-multiscen "12:34:56"` (the string argument is the estimated time needed for each scenario)
 * Analyze outputs with `simul/go-hpc-analyze "12:34:56"`
 * Merge the results into one RData file for future plots: `Rscript multiscen-merge-results.R`
 * Plot figures: `Rscript figures-MAIN.R` and `Rscript figures-SI.R`




