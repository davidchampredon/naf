# naf
## New agent-based model for flu

Epidemiological model is coded in C++ and wrapped in R for output analysis. 

Preferred way to run the simulations:
* 0 [Optional]) Make sure the C++ code is correct by compiling it: `make naf`. Code can be run with `./naf`.
* 1) Build the R `naf` library: in folder `Rlibrary`, execute `./build_library`. If there is a problem at this stage only (if the previous optional stage is run without any error), the problem probably lies in `Rwrap.cpp`.
* 2 [Optional]) Fit the age distributions conditional on household size on the global age distribution: in folder `fit`, run `go-fit-hhszage`. Once the fit is done, generate the best fitted distributions by running the script `go-gen-ad-hhsz` in the folder `data/households`. The best fitted parameters must be copy-pasted in `gen-ad-hhsz.R` (TO DO: change this).
* 3) Simulations with scenario comparison are launched with the script `go-simul` in the `simul` folder.  



