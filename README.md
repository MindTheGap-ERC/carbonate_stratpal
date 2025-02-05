# lineage_ecol
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


## Repository structure  

* _code_ : folder with R code
  * _make_table_and_plots.R_ : script to generate table 1 and figures
  * _make_maps.R_ : script to generate the maps
  * _download_data.R_ : script to download data from Zenodo
* _data_ : folder for raw data. Initially empty, will be filled with downloaded data after the script in `code/download_data.R` is run.
* _figs_ : folder for figures. Initially empty
* _renv_ : folder used by the `renv` package
* _.gitignore_ : untracked files
* _.Rprofile_ : session info
* _messinian_db.Rproj_ : RProject file
* _README_ : README file
* _renv.lock_ : lock file for `renv` package

## References

This repository uses data from
* site the CarboKitten webpage

It uses sea level data from 
* Holland, S. M., & Patzkowsky, M. E. (2015). The stratigraphy of mass extinction. Palaeontology, 58(5), 903-924.


## Funding

This work was supported by Utrecht University.
