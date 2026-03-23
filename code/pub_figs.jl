using CarboKitten.Visualization
import CarboKitten.Visualization: summary_plot
using CarboKitten.Export: read_header, read_volume, read_slice, group_datasets
using CarboKitten.Utility: in_units_of
using HDF5
using Unitful
using CairoMakie
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution

const na = [CartesianIndex()]

elevation(h, d) =
    let bl = h.initial_topography[d.slice..., na],
        sr = h.axes.t[end] * h.subsidence_rate

        bl .+ d.sediment_thickness .- sr
    end

water_depth(header, data) =
    let h = elevation(header, data),
        wi = data.write_interval,
        s = header.subsidence_rate .* (header.axes.t[1:wi:end] .- header.axes.t[end]),
        l = header.sea_level[1:wi:end]

        h .- (s.+l)[na, :]
    end

function dominant_facies!(ax, header, data;
    smooth_size::NTuple{2,Int}=(3, 11),
    colors=Makie.wong_colors())
    n_facies = size(data.production)[1]
    colormax(d) = getindex.(argmax(d; dims=1)[1, :, :], 1)

    wi = data.write_interval
    dominant_facies = colormax(data.deposition)
    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
    wd = zeros(Float64, length(header.axes.x), length(header.axes.t[1:wi:end]))
    blur(water_depth(header, data) / u"m", wd)

    ax.ylabel = "time [Myr]"
    ax.xlabel = "position [km]"

    xkm = header.axes.x |> in_units_of(u"km")
    tmyr = header.axes.t[1:wi:end] |> in_units_of(u"Myr")
    ft = heatmap!(ax, xkm, tmyr, dominant_facies;
        colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
        colorrange=(0.5, n_facies + 0.5))
    contourf!(ax, xkm, tmyr, wd;
        levels=[0.0, 10000.0], colormap=Reverse(:grays))
    #contour!(ax, xkm, tmyr, wd;
    #    levels=[0], color=:black, linewidth=2)
    return ft
end

fig = Figure(size=(1200, 1000), backgroundcolor=:gray80)

TAG1 = "platform"
TAG2 = "ramp"

header, data = read_slice("data/$(TAG1).h5", :profile)
ax = Axis(fig[1:2, 1:2])
sediment_profile!(ax, header, data, show_unconformities = true, show_coeval_lines = false)
#coeval_lines!(ax, header, data, [0.25:0.5:3.75;]u"Myr", linewidth = 3, color = :black, linestyle = :solid)

header, data = read_slice("data/$(TAG2).h5", :profile)
ax = Axis(fig[1:2, 4:5])
sediment_profile!(ax, header, data, show_unconformities = true, show_coeval_lines = false)
#coeval_lines!(ax, header, data, [0.25:0.5:3.75;]u"Myr", linewidth = 3, color = :black, linestyle = :solid)

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

ax = Axis(fig[3, 4:5])
header, data = read_slice("data/$(TAG2).h5", :profile)
dominant_facies!(ax, header, data)

ax = Axis(fig[3, 1:2])
header, data = read_slice("data/$(TAG1).h5", :profile)
dominant_facies!(ax, header, data)

ax = Axis(fig[3, 3], title = "sea level curve",xlabel = "sealevel [m]",
limits = (nothing, (header.axes.t[1] |> in_units_of(u"Myr"), header.axes.t[end]|> in_units_of(u"Myr"))))
lines!(ax, header.sea_level |> in_units_of(u"m") , header.axes.t|> in_units_of(u"Myr"))

save("figs/ms/fig1.png", fig)
