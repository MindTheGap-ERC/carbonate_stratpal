# Diversity record of fossils with different ecological niches

BSc thesis of Anna Jansen

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
instantiate
```

this will download and precompile all dependencies. It will take a while when you run it the first time, but each next time you run a Julia script in this repo will be much faster.

4. Exit the package mode by pressing the backspace button. This will bring you back to the normal Julia prompt.

5. Run a script using the command `include`, for example:

```{julia}
include("src/Run_model.jl")
```


### Repository structure
*Folder structure goes here

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

## Copyright

Copyright 2025 Anna Jansen and Utrecht University

## License

Apache 2.0 License, see LICENSE file for license text.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.
![European Union and European Research Council logos](https://erc.europa.eu/sites/default/files/2023-06/LOGO_ERC-FLAG_FP.png)
