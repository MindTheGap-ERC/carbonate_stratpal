# ~/~ begin <<docs/src/model-alcap.md#examples/model/alcap/plot.jl>>[init]
#| creates: docs/src/_fig/alcaps-alternative.png
#| requires: data/output/alcap-example.h5
#| collect: figures
import Pkg
Pkg.add("GLMakie")

using GLMakie
using CarboKitten.Visualization

save("figs/lineage_ecology.png", summary_plot("data/alcap-example.h5"))
# ~/~ end
