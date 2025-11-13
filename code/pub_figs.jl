using CarboKitten.Visualization
import CarboKitten.Visualization: summary_plot
using CarboKitten.Export: read_header, read_volume, read_slice, group_datasets
using CarboKitten.Utility: in_units_of
using HDF5
using Unitful
using CairoMakie

fig = Figure(size=(1200, 1000), backgroundcolor=:gray80)

TAG1 = "platform"
TAG2 = "ramp"

header, data = read_slice("data/$(TAG1).h5", :profile)
ax = Axis(fig[1:2, 1:2])
sediment_profile!(ax, header, data, show_unconformities = true, show_coeval_lines = false)
coeval_lines!(ax, header, data, [0.25:0.5:3.75;]u"Myr", linewidth = 3, color = :black, linestyle = :solid)

header, data = read_slice("data/$(TAG2).h5", :profile)
ax = Axis(fig[1:2, 4:5])
sediment_profile!(ax, header, data, show_unconformities = true, show_coeval_lines = false)
coeval_lines!(ax, header, data, [0.25:0.5:3.75;]u"Myr", linewidth = 3, color = :black, linestyle = :solid)

ax = Axis(fig[1, 3])
max_depth = minimum(header.initial_topography)
h5open("data/$(TAG1).h5", "r") do fid
    production_curve!(ax, fid["input"], max_depth=max_depth)
end

ax = Axis(fig[2, 3])
max_depth = minimum(header.initial_topography)
h5open("data/$(TAG2).h5", "r") do fid
    production_curve!(ax, fid["input"], max_depth=max_depth)
end

ax1 = Axis(fig[3, 1])
ax2 = Axis(fig[3, 2])
header, data = read_slice("data/$(TAG1).h5", :profile)
wheeler_diagram!(ax1, ax2, header, data)

ax1 = Axis(fig[3, 4])
ax2 = Axis(fig[3, 5])
header, data = read_slice("data/$(TAG2).h5", :profile)
wheeler_diagram!(ax1, ax2, header, data)

ax = Axis(fig[3, 3], title = "sea level curve",xlabel = "sealevel [m]",
limits = (nothing, (header.axes.t[1] |> in_units_of(u"Myr"), header.axes.t[end]|> in_units_of(u"Myr"))))
lines!(ax, header.sea_level |> in_units_of(u"m") , header.axes.t|> in_units_of(u"Myr"))

save("figs/test.png", fig)
