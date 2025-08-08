# lineage_ecol

The supplementary material for "title of research". All code is included to produce this research's figures, plots and statistics. Not all produced figures are presented in the final paper. Additional figures are included in this repository.

## Authors

__Sidney Bickerton__  (maintainer, contributor)  
Utrecht University  
email: s.j.bickerton [at] uu.nl  

__Niklas Hohmann__  (maintainer, contributor)  
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

## Purpose
The purpose of this research is to employ forward modelling to predict the location of stratigraphically created last occurrence clusters and distinguish them from clusters created by increased extinction rate on a carbonate platform. This research also aims to determine differences in patterns of last occurrence clusters in carbonate platforms with patterns found in siliciclastic platforms.    

## License

Copyright 2025 Sidney Bickerton, Utrecht University and the Netherlands eScience Center

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

## Requirements

Base R (version >= 4.0) and the RStudio IDE.

## Reproduction

To reproduce the data, first use `code/run.jl` in Julia. This creates the .csv files `lineage-example_adm`, `lineage-example_sac`, and `lineage-example_sc`, containing the data used for all the scenarios. It also creates the file `lineage-example.h5`. This file is used in `code/plot.jl` in Julia to produce the figure of the carbonate platform which is saved as `lineage_ecology.png`. 
Then run `code/package_installation.R` in RStudio to install `admtools` and `StratPal`. This makes you able to simulate Age Depth Models and the phenotypic evolution behind them. Then run `code/make_figures.R` which uses the raw data from `CarboKitten.jl` to produce all the figures used in this research. 

## Repository structure  

* _code_ : folder with R code
  * _run.jl_ : creates the raw data in a `.csv` file for the carbonate platform using `CarboKitten.jl`
  * _plot.jl_ : creates a plot of the carbonate platform using `CarboKitten.jl`
  * _make_figures.R_ : script to generate all figures
  * _package_installation.R_ : script to install the used packages and download data
  * _stats_on_data.R_ : script to calculate statistics and specific numbers mentioned in the paper
* _data_ : folder with raw data from CarboKitten.
* _figs_ : folder for all used figures figures. 
* _.gitignore_ : untracked files.
* _LICENSE_ : licence of research.
* _README_ : README file.

## References

This repository uses data from
* the [`CarboKitten.jl`](https://mindthegap-erc.github.io/CarboKitten.jl/) webpage

Carbonate production data is taken from
* Burgess, P. M. (2013). CarboCAT: A cellular automata model of heterogeneous carbonate strata. Computers & geosciences, 53, 129-140.
Eustatic sea level data it taken from 
* Holland, S. M., & Patzkowsky, M. E. (2015). The stratigraphy of mass extinction. Palaeontology, 58(5), 903-924.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.

