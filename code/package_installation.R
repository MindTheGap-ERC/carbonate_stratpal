install.packages("remotes") # can install from github
install.packages("paleoTS") # paleontological time series

remotes::install_github(repo = "MindTheGap-ERC/admtools",
                        build_vignettes = TRUE,
                        dependencies = TRUE)

remotes::install_github(repo = "MindTheGap-ERC/StratPal",
                        build_vignettes = TRUE,
                        ref = "dev",
                        dependencies = TRUE)

library(StratPal) # biology
library(admtools) # age-depth
library(paleoTS)  # time series analysis

#Read data from CarboKitten
data_kitten = read.csv("C:/Users/sidne/OneDrive/Documenten/AA Utrecht/Guided research/Git repositories/Stratigraphic Paleobiology for Phenotypic Evolution/CarboKitten.jl/data/output/alcap-example3_adm.csv")