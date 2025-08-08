import Pkg
Pkg.add("GLMakie")

using GLMakie
using CarboKitten.Visualization

save("figs/lineage_ecology.png", summary_plot("data/lineage-example.h5"))
# ~/~ end
