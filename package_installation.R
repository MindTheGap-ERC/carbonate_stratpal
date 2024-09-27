remotes::install_github(repo = "MindTheGap-ERC/admtools",
                        build_vignettes = TRUE,
                        ref = "dev",
                        dependencies = TRUE)
remotes::install_github(repo = "MindTheGap-ERC/StratPal",
                        build_vignettes = TRUE,
                        ref = "dev",
                        dependencies = TRUE)

library(StratPal)
library(admtools)

adm = tp_to_adm(scenarioA$t_myr, scenarioA$h[,"2km"])

plot(adm)

x = stasis_sl(scenarioA$t_myr)

plot(reduce_to_paleoTS.pre_paleoTS(x))

a = x |> time_to_strat(adm) |> reduce_to_paleoTS.pre_paleoTS() |> plot(xlab = "Stratigraphic Position")


