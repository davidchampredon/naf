# naf
## New agent-based model for flu

Epidemiological model is coded in C++ and wrapped in R for output analysis. 

Preferred way to run the simulations:
* 1) Make sure the C++ code is correct by compiling it: `make naf`. Code can be run with `./naf`.
* 2) Build the R `naf` library: in folder `Rlibrary`, execute `./build_library`. If there is a problem at this stage only, the problem probably lies in `Rwrap.cpp`.
* 3) The R script `test_naf_run.R` defines inputs, processes and plots outputs of the simulation. Execute `Rscript test_naf_run.R` and a file `plot_TEST_naf.pdf` is created plpotting a summary of the simulation.



