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

## License



## Requirements

Base R (version >= 4.0) and the RStudio IDE.

## Reproduction

To reproduce the data, first use run.jl in julia. Then run package_installation.R in Rstudio to install admtools and StratPal which makes you able to simulate age depth models and the phenotypic evolution behind them. Then run make_figures which uses the raw data from CarboKitten to produce all the figures used in this research. 

## Repository structure  

* _code_ : folder with R code
  * _make_figures.R_ : script to generate all figures
  * _package_installation.R_ : script to install the used packages and download data
  * _stats_on_data.R_ : script to calculate statistics and additional figures
  * _run.jl_ : creates the raw data in a .csv file for the carbonate platform using CarboKitten
  * _plot.jl_ : creates a plot of the carbonate platform using CarboKitten
* _data_ : folder with raw data from CarboKitten.
* _figs_ : folder for all used figures figures. 
* _.gitignore_ : untracked files.
* _LICENSE_ : licence of research.
* _README_ : README file.


## References

This repository uses data from
* site the CarboKitten webpage

It uses sea level data from 
* Holland, S. M., & Patzkowsky, M. E. (2015). The stratigraphy of mass extinction. Palaeontology, 58(5), 903-924.


## Funding

This work was supported by Utrecht University.
