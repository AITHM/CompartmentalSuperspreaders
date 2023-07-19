# CompartmentalSuperspreaders

This repository contains all of the code and data required to generate the outputs in the article titled: "Replicating superspreader dynamics with compartmental models".

The main analysis has been written in the [Julia](https://julialang.org) (v1.8) programming language, whilst post-processing and plotting is performed in [R](https://www.r-project.org) (v4.2.3).

## Dependencies
The explicit dependencies for the Julia analysis can be found in `Manifest.toml`, whilst the post-processing in R uses

	- tidyverse v2.0.0
	- ggsci v3.0.0
	- xtable v1.8.4

## Repository contents
	- data/offspring: offspring distribution data used for model fitting
	- outputs: post-processing utilities and empty directories to store computational outputs (including both Julia and R generated files)
	- params: baseline and sensitivity parameters
	- src: Julia source code for model fitting and analysis
	- `run_offspring.jl`: main analysis script to reproduce article results

## Reproducing results
### Simulation
To reproduce the computational results of the article navigate to the repository directory in Julia using the `cd("...filepath...")` command:
```julia
julia> cd("...filepath...")
```

From there, activate the project environment by first entering `]` to enter package mode, and then activating all:
```julia
julia> ]
pkg> activate .
(CompartmentalSuperspreaders) pkg> instantiate
```
To escape package mode and return to the normal prompt simply press backspace.

Finally, to run the main analysis (along with all sensitivity simulations) include the `run_offspring.jl` file:
```julia
julia> include("run_offspring.jl")
```

Running this script will populate `outputs\offspring` and its sub-directories.


### Plotting
To plot the outputs navigate to the repository directory in RStudio using the `setwd(...filepath...)` command, and then open the `outputs\offspring` directory in the script editor.

To generate all baseline plots and tables simply select all (Ctrl + A) and run (Ctrl + Enter). Sensitivity analysis plots can be created by changing the definition of the `output_dir` variable given on line 7. In all cases, the figures are saved as .png images in the corresponding `outputs/offspring` directory.
Note that for windows, the `plot_results.R` file cannot be run via source due to the presence of unicode characters and this may also influence the appearance of some plots. 



