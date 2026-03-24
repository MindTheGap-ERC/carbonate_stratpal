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

TAG1 = "platform"
TAG2 = "ramp"

fig = Figure(size=(1200, 1000), backgroundcolor=:gray80)

ax = Axis(fig[1:2, 1:2])
header, data = read_slice("data/$(TAG1).h5", :profile)
dominant_facies!(ax, header, data)

st_bdry = 0.25:0.5:3.75
text = ["TST1", "HST1", "RST1", "LST1", "TST2", "HST2", "RST2", "LST2", "TST3"]
sampling_locations = [3,7.5,10.5, 12, 18]

x = header.axes.x |> in_units_of(u"km")
x_max = x[end]
x_min = x[1]
t = header.axes.t |> in_units_of(u"Myr")
t_max = t[end]
t_min = t[1]
st_bdr_ext = [t_min, st_bdry..., t_max]
plot_pos = 0.5 .* (st_bdr_ext[1:end-1] + st_bdr_ext[2:end])


for pos in st_bdry
    lines!(ax, [x_min, x_max], repeat([pos], 2), color = :black, linewidth = 5)
end


for pos in sampling_locations
    lines!(ax, [pos,pos],[0,4], linestyle = :dash, linewidth = 5)
end
for (i, pos) in enumerate(plot_pos)
    text!(ax, 5, pos, text = text[i])
end



save("figs/sm/sfig3_wheeler_$TAG1.png", fig)


ax = Axis(fig[1:2, 1:2])
header, data = read_slice("data/$(TAG2).h5", :profile)
dominant_facies!(ax, header, data)

st_bdry = 0.25:0.5:3.75
text = ["TST1", "HST1", "RST1", "LST1", "TST2", "HST2", "RST2", "LST2", "TST3"]
sampling_locations = [3,7.5,10.5, 12, 18]

x = header.axes.x |> in_units_of(u"km")
x_max = x[end]
x_min = x[1]
t = header.axes.t |> in_units_of(u"Myr")
t_max = t[end]
t_min = t[1]
st_bdr_ext = [t_min, st_bdry..., t_max]
plot_pos = 0.5 .* (st_bdr_ext[1:end-1] + st_bdr_ext[2:end])


for pos in st_bdry
    lines!(ax, [x_min, x_max], repeat([pos], 2), color = :black, linewidth = 5)
end


for pos in sampling_locations
    lines!(ax, [pos,pos],[0,4], linestyle = :dash, linewidth = 5)
end
for (i, pos) in enumerate(plot_pos)
    text!(ax, 5, pos, text = text[i])
end



save("figs/sm/sfig4_wheeler_$TAG2.png", fig)