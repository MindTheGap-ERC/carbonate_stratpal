module Plot_data

using CarboKitten
using Unitful
using GLMakie
using CarboKitten.Visualization: summary_plot

const TAG = "example_Anna"
const PATH = "data"

function plot_model()
    summary_plot("data/example_Anna.h5")
end

end

Plot_data.plot_model()