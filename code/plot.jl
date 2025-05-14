#| creates: docs/src/_fig/alcaps-alternative.png
#| requires: data/output/lineage-example.h5

import Pkg
Pkg.add("GLMakie")

using GLMakie
using CarboKitten.Visualization

save("figs/lineage_ecology.png", summary_plot("data/lineage-example.h5"))
# ~/~ end
