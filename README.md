# Stratigraphic paleobiology of carbonate systems

Build upon the work of Anna Jansen and Sidney Bickerton

## Authors

__Sidney Bickerton__
Utrecht University  
email: s.j.bickerton [at] uu.nl  

__Anna Jansen__

__Niklas Hohmann__  
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

__Emilia Jarochowska__  
Utrecht University  
email: e.b.jarochowska [at] uu.nl  
Web page: [www.uu.nl/staff/EBJarochowska](https://www.uu.nl/staff/EBJarochowska)  
ORCID: [0000-0001-8937-9405](https://orcid.org/0000-0001-8937-9405)

__Xianyi Liu__  
Utrecht University  
email: x.liu6 [at] uu.nl  
Web page: [www.uu.nl/staff/XLiu6](https://www.uu.nl/staff/XLiu6)  
ORCID: [0000-0002-3851-116X](https://orcid.org/0000-0002-3851-116X)


## Requirements

Base R (version >= 4.0) and the RStudio IDE, Julia >= 1.10



### Installation and running instructions

#### Running a CarboKitten model

You need to have Julia. Follow [the instructions on its webpage](https://julialang.org/downloads/).

1. Clone this repository and go to its folder

2. Open Julia (either in your terminal or in VS Code)

3. Initialize all dependencies by typing the following:

```{julia}
using Pkg
```

press `]` to enter the package mode and there type:

```{julia}
activate .
instantiate
```

this will download and precompile all dependencies. It will take a while when you run it the first time, but each next time you run a Julia script in this repo will be much faster.

4. Exit the package mode by pressing the backspace button. This will bring you back to the normal Julia prompt.

5. Run a script using the command `include`, for example:

```{julia}
include("src/run.jl")
```
to generate the data,

```{julia}
include("src/plot.jl")
```

to generate the transect plot, and 

```{julia}
include("src/extract_data.jl")
```

to extract water depth & age-depth models as .csv files.

After generating the data, you can run the analysis in R using

```{R}
source("code/analysis.R")
```

This will generate figures under `figs/`

## Repository structure

* code: Julia & R scripts
* data: model outputs. Initially empty, filled after CarboKitten is run.
* figs: figures. Initially empty, filled after CarboKitten/the R code is run
* .gitignore: untracked files
* carbonate_stratpal.Rproj: R project file
* CONTRIBUTING.md: contribution guidelines
* LICENSE: Apache 2.0 License text
* Manifest.toml: Julia project infrastructure
* Project.toml: Julia project infrastructure
* README.md: this file

## Copyright

Copyright 2025 Anna Jansen, Sidney Bickerton, Utrecht University and the Netherlands eScience Center

## License

Apache 2.0 License, see LICENSE file for license text.

## References

This repository uses data from
* the [`CarboKitten.jl`](https://mindthegap-erc.github.io/CarboKitten.jl/) webpage

Carbonate production data is taken from
* Burgess, P. M. (2013). CarboCAT: A cellular automata model of heterogeneous carbonate strata. Computers & geosciences, 53, 129-140.
Eustatic sea level data it taken from 
* Holland, S. M., & Patzkowsky, M. E. (2015). The stratigraphy of mass extinction. Palaeontology, 58(5), 903-924.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.
![European Union and European Research Council logos](https://erc.europa.eu/sites/default/files/2023-06/LOGO_ERC-FLAG_FP.png)

