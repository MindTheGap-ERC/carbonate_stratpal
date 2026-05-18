# Stratigraphic paleobiology of carbonate systems

Build upon the work of Anna Jansen and Sidney Bickerton

## Authors

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

__Sidney Bickerton__

__Anna Jansen__

## Usage

For reproducibility instructions, see REPRODUCEME.md

## Requirements

Base R v4.5.3 and the RStudio IDE, Julia 1.12.6

## Repository structure

* code: Julia & R scripts
* data: model outputs. Initially empty, filled after CarboKitten is run.
* figs: figures. Initially empty, filled after CarboKitten/the R code is run
    * ms/: figures for manuscript
    * sm/: figures for supplementary materials
* renv: `renv` setup
* .gitignore: untracked files
* .Rprofile: R session/project setup
* carbonate_stratpal.Rproj: R project file
* CONTRIBUTING.md: contribution guidelines
* LICENSE: Apache 2.0 License text
* Manifest.toml: Julia project infrastructure
* Project.toml: Julia project infrastructure
* README.md: you are here
* renv.lock: lock file for `renv`
* REPRODUCEME.md: Instructions for reproduction

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

