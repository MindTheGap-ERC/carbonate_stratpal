# Instructions for reproduction

## Simulating Carbonate Platforms

**Required data**: None
**Generated data**: Carbonate platform data under `data/` as `.h5` or `.csv` files, and summary figures under `figs/`

After installing Julia via [Juliaup](https://github.com/JuliaLang/juliaup), run

```julia
julia +1.12.6 --project=. -i -e "using Pkg; Pkg.instantiate()"
```

in the command line at the root of this repository. This will start a Julia REPL session, activate the project, and install and precompile all dependencies.

Then run

```Julia
include("code/run_all.jl")
```

to run all carbonate platform simulations, extract data as .csv from them, and plot the figures. You can also run and inspect individual components by including them individually.

## Stratigraphic Paleobiology Simulations

**Required data:** Carbonate platform data as generated above
**Generated data:** Figures under `figs/`

Open the file `carbonate_stratpal.Rproj` in the RStudio IDE. This will open the project, set all paths correctly, and install the `renv` package (if not already installed). Next, run

```R
renv::restore()
```

to install all required packages in their correct versions.

Runing

```R
source("code/analysis.R")
```

in the console will produce all figures under `figs/`.