install.packages("remotes") # can install from github
install.packages("paleoTS") # paleontological time series



remotes::install_github(repo = "MindTheGap-ERC/admtools",
                        build_vignettes = TRUE,
                        dependencies = TRUE)

remotes::install_github(repo = "MindTheGap-ERC/StratPal",
                        build_vignettes = TRUE,
                        ref = "dev",
                        dependencies = TRUE)

library(StratPal) # biology part
library(admtools) # age-depth
library(paleoTS)  # time series analysis

x = paleoTS::sim.GRW()
plot(x)

# fit models
fit3models(x)
fit9models(x)

?`paleoTS-package`
fitSimple(x, model = "OU")
# fits only one model, have to do model comparison afterwards
# inspect the code of fit3models
fit3models # ctl + enter


## new stratpal functionality

adm = tp_to_adm(scenarioA$t_myr, scenarioA$h[,"2km"])

plot(adm)
# paleots is a summary

# pre-paleoTS

x = stasis_sl(1:3) # sl stand for specimen level
x = random_walk_sl()
?random_walk_sl

x$t # times
x$vals # individual trait values

# plotting pipeline
seq(0.1, 1.9, by  = 0.1) |>
  random_walk_sl() |>
  time_to_strat(adm) |>
  reduce_to_paleoTS() |>
  plot()

# use `reduce_to_paleoTS` before plotting
?reduce_to_paleoTS

# analysis pipeline
seq(0.1, 1.9, by  = 0.1) |>
  random_walk_sl() |>
  time_to_strat(adm) |>
  reduce_to_paleoTS() |>
  fit3models()
test test
