module Plot
import CairoMakie
using CarboKitten.Algorithms: skeleton
import CarboKitten.Visualization: summary_plot, sediment_profile, sediment_profile!, production_curve!
using CarboKitten.Export: Header, Data, DataSlice, read_header, read_volume, read_slice, group_datasets
using CarboKitten.Utility: in_units_of
using CarboKitten.Output.Abstract: stratigraphic_column
using HDF5
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution
using CairoMakie
import CairoMakie


const na = [CartesianIndex()]

const sampl_loc = [3, 7.5, 10.5, 12, 18]u"km"
const st_div = 0.25u"Myr":0.5u"Myr":3.75u"Myr"
const st_labels = ["TST1", "HST1", "RST1", "LST1", "TST2", "HST2", "RST2", "LST2", "TST3"]
const st_pos = [0.125u"Myr", 0.5u"Myr":0.5u"Myr":3.5u"Myr"..., 3.875u"Myr"]
const pos_linestyle=:solid
const pos_linewidth=2
const pos_color=:black
const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")

const time_axis_label = "Elapsed model time [Myr]"
const water_depth_label = "Water depth [m]"
const production_label = "Production [m/Myr]"
const distance_label = "Distance from shore [km]"
const depth_label = "Depth [m]"
const subfig_label_fontsize = 24

TAG1 = "platform"
TAG2 = "ramp"

function explode_quad_vertices(v::Array{Float64,3})
    w, h, d = size(v)
    points = zeros(Float64, w, h - 1, 2, d)
    n_vertices = 2 * w * (h - 1)
    n_quads = (w - 1) * (h - 1)
    @views points[:, :, 1, :] = v[1:end, 1:end-1, :]
    @views points[:, :, 2, :] = v[1:end, 2:end, :]
    idx = reshape(1:n_vertices, w, (h - 1), 2)
    vtx1 = reshape(idx[1:end-1, :, 1], n_quads)
    vtx2 = reshape(idx[2:end, :, 1], n_quads)
    vtx3 = reshape(idx[2:end, :, 2], n_quads)
    vtx4 = reshape(idx[1:end-1, :, 2], n_quads)
    return reshape(points, n_vertices, d),
           vcat(hcat(vtx1, vtx2, vtx3), hcat(vtx1, vtx3, vtx4))
end

function profile_plot!(f::F, ax::Axis, header::Header, data::DataSlice; mesh_args...) where {F}
    color = f.(eachslice(data.deposition, dims=(2, 3)))
    profile_plot!(ax, header, data; color=color, mesh_args...)
end

function plot_unconformities(ax::Axis, header::Header, data::DataSlice, minwidth::Int; kwargs...)
    n_facies, n_x, n_t = size(data.production)
    total_subsidence = (header.axes.t[end] - header.axes.t[1]) * header.subsidence_rate
    initial_topography = header.initial_topography[data.slice...]
    sc = stratigraphic_column(data)
    h = repeat(initial_topography .- total_subsidence, 1, n_t+1)
    @views h[:, 2:end] .+= cumsum(sum(sc, dims=1)[1,:,:], dims=2)
    x = header.axes.x |> in_units_of(u"km")
    combined_acc = dropdims(sum(sc, dims = 1), dims = 1) |> in_units_of(u"m")
    wi = data.write_interval
    hiatus = skeleton(combined_acc .< 0.000001, minwidth=minwidth)

    if !isempty(hiatus[1])
        verts = [(x[pt[1]], h[pt...] |> in_units_of(u"m")) for pt in hiatus[1]]
        linesegments!(ax, vec(permutedims(verts[hiatus[2]])); kwargs...)
    end
end

function profile_plot!(ax::Axis, header::Header, data::DataSlice; color::AbstractArray, mesh_args...)
    x = header.axes.x |> in_units_of(u"km")
    t = header.axes.t |> in_units_of(u"Myr")

    n_facies, n_x, n_t = size(data.production)
    total_subsidence = (header.axes.t[end] - header.axes.t[1]) * header.subsidence_rate
    initial_topography = header.initial_topography[data.slice...]
    sc = stratigraphic_column(data)
    h = repeat(initial_topography .- total_subsidence, 1, n_t+1)
    @views h[:, 2:end] .+= cumsum(sum(sc, dims=1)[1,:,:], dims=2)

    verts = zeros(Float64, n_x, n_t+1, 2)
    @views verts[:, :, 1] .= x
    @views verts[:, :, 2] .= h |> in_units_of(u"m")
    v, f = explode_quad_vertices(verts)

    total_subsidence = header.subsidence_rate * (header.axes.t[end] - header.axes.t[1])
    bedrock = (header.initial_topography[data.slice...] .- total_subsidence) |> in_units_of(u"m")
    lower_limit = minimum(bedrock) - 20
    band!(ax, x, lower_limit, bedrock; color=:gray, label="initial topography")
    lines!(ax, x, bedrock; color=:black, label="initial topography")
    ylims!(ax, lower_limit + 10, nothing)
    xlims!(ax, x[1], x[end])
    ax.xlabel = "Distance from shore [km]"
    ax.ylabel = "Depth [m]"

    c = reshape(color, n_x * n_t)
    mesh!(ax, v, f; color=vcat(c, c), mesh_args...)
end

function plot_sealevel!(ax::Axis, header::Header)
    sea_level = header.sea_level[end] |> in_units_of(u"m")
    hlines!(ax, sea_level, color=:lightblue, linewidth=5, label="end sea level")
end



function sediment_profile!(ax::Axis, header::Header, data::DataSlice;
                           show_unconformities::Union{Nothing,Bool,Int} = true,
                           show_coeval_lines::Union{Bool,Tuple{Int, Int},Vector{Int},Vector{Time}} = true,
                           show_sealevel::Bool = true,
                           title = "Sediment profile")
    x = header.axes.x |> in_units_of(u"km")
    t = header.axes.t |> in_units_of(u"Myr")

    n_facies, n_x, n_t = size(data.production)
    total_subsidence = (header.axes.t[end] - header.axes.t[1]) * header.subsidence_rate
    initial_topography = header.initial_topography[data.slice...]
    sc = stratigraphic_column(data)
    h = repeat(initial_topography .- total_subsidence, 1, n_t+1)
    @views h[:, 2:end] .+= cumsum(sum(sc, dims=1)[1,:,:], dims=2)

    if show_sealevel
        plot_sealevel!(ax, header)
    end

    plot = profile_plot!(argmax, ax, header, data; alpha=1.0,
        colormap=cgrad(Makie.wong_colors()[1:n_facies], n_facies, categorical=true))

    minwidth = 10
    if show_unconformities
        plot_unconformities(ax, header, data,minwidth; label = "unconformities",
                        color=:white, linestyle=:dash, linewidth=1)
    end

    if show_coeval_lines
        coeval_lines!(ax, header, data, age_depth_model(header, data), collect(st_div); color=:black, linestyle=:dash)
    end

    ax.title = title
    return plot
end

elevation(h::Header, d::DataSlice) =
    let bl = h.initial_topography[d.slice..., na],
        sr = h.axes.t[end] * h.subsidence_rate

        bl .+ d.sediment_thickness .- sr
    end

water_depth(header::Header, data::DataSlice) =
    let h = elevation(header, data),
        wi = data.write_interval,
        s = header.subsidence_rate .* (header.axes.t[1:wi:end] .- header.axes.t[end]),
        l = header.sea_level[1:wi:end]

        h .- (s.+l)[na, :]
    end


function dominant_facies!(ax::Axis, header::Header, data::DataSlice;
    show::Symbol = :both,
    smooth_size::NTuple{2,Int} = (3, 11),
    colors = Makie.wong_colors())

    if show ∉ [:model, :preserved, :both]
        error("expected argument `show` to be one of `:model`, `:preserved`, `:both`; got $(show)")
    end

    prec = 10^-8 # prec. in m - below acc. is considered 0

    n_facies = size(data.production)[1]
    colormax(d) = getindex.(argmax(d; dims=1)[1, :, :], 1)

    wi = data.write_interval
    blur = convolution(Shelf, ones(Float64, smooth_size...) ./ *(smooth_size...))
    wd = zeros(Float64, length(header.axes.x), length(header.axes.t[1:wi:end]))
    blur(water_depth(header, data) / u"m", wd)

    ax.ylabel = "Elapsed model time [Myr]"
    ax.xlabel = "Distance from shore [km]"
    ax.title = "Wheeler diagram"
    
    xkm = header.axes.x |> in_units_of(u"km")
    tmyr = header.axes.t[1:wi:end] |> in_units_of(u"Myr")


    ft = if show == :model
        dominant_facies = colormax(data.deposition)
        dominant_facies = Matrix{Union{Missing, Int}}(dominant_facies)
        dominant_facies[ wd .> 0] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies;
            colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
            colorrange=(0.5, n_facies + 0.5),
            nan_color = :white)
    elseif show == :preserved
        sc = stratigraphic_column(data)
        dominant_facies = colormax(sc)
        dominant_facies = Matrix{Union{Missing, Int}}(dominant_facies)
        combined_acc = dropdims(sum(sc, dims = 1), dims = 1) |> in_units_of(u"m")
        dominant_facies[ combined_acc .< prec] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies;
            colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
            colorrange=(0.5, n_facies + 0.5),
            nan_color = :white)
    else
        sc = stratigraphic_column(data)
        dominant_facies_model = colormax(data.deposition)
        dominant_facies_model = Matrix{Union{Missing, Int}}(dominant_facies_model)
        dominant_facies_model[ wd .> 0] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies_model;
            colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
            colorrange=(0.5, n_facies + 0.5), alpha=0.3,
            nan_color = :white)
        dominant_facies_preserved = colormax(sc)
        dominant_facies_preserved = Matrix{Union{Missing, Int}}(dominant_facies_preserved)
        combined_acc = dropdims(sum(sc, dims = 1), dims = 1) |> in_units_of(u"m")
        dominant_facies_preserved[ combined_acc .< prec] .= missing
        heatmap!(ax, xkm, tmyr, dominant_facies_preserved;
            colormap=cgrad(colors[1:n_facies], n_facies, categorical=true),
            colorrange=(0.5, n_facies + 0.5),
            nan_color = :transparent)
    end

    #contourf!(ax, xkm, tmyr, wd;
    #    levels=[0.0, 10000.0], colormap=Reverse(:grays))
    #contour!(ax, xkm, tmyr, wd;
    #    levels=[0], color=:black, linewidth=2)
    return ft
end

function coeval_lines!(ax::Axis, header::Header, data::DataSlice, adm::AbstractMatrix{Amount}, tics::Vector{Int}; kwargs...)
    x = header.axes.x |> in_units_of(u"km")
    h = adm |> in_units_of(u"m")
    for t in tics
        lines!(ax, x, h[:, t]; kwargs...)
    end
end

function coeval_lines!(ax::Axis, header::Header, data::DataSlice, adm::AbstractMatrix{Amount}, tics::Vector{Time}; kwargs...)
    t = header.axes.t[1:data.write_interval:end]
    indices::Vector{Int} = [searchsortedfirst(t, i) for i in tics]
    coeval_lines!(ax, header, data, adm, indices; kwargs...)
end

function age_depth_model(header::Header, data::DataSlice)
    x = header.axes.x |> in_units_of(u"km")
    t = header.axes.t |> in_units_of(u"Myr")

    n_facies, n_x, n_t = size(data.production)
    total_subsidence = (header.axes.t[end] - header.axes.t[1]) * header.subsidence_rate
    initial_topography = header.initial_topography[data.slice...]
    sc = stratigraphic_column(data)
    h = repeat(initial_topography .- total_subsidence, 1, n_t+1)
    @views h[:, 2:end] .+= cumsum(sum(sc, dims=1)[1,:,:], dims=2)
    return h
end

plot_supp_fig = function(tag, name)
    fig = Figure(size=(1200, 1000), backgroundcolor=:gray80)

    header, data = read_slice("data/$(tag).h5", :profile)

    adm = age_depth_model(header, data)


    ax_sed = Axis(fig[1,1])
    ax_wheel = Axis(fig[2,1])
    linkxaxes!(ax_sed, ax_wheel)
    sediment_profile!(ax_sed, header, data, show_coeval_lines=false, show_unconformities = true)
    coeval_lines!(ax_sed, header, data, adm, collect(st_div); color=:black, linestyle=:dash)
    for loc in sampl_loc
        lines!(ax_sed,fill(loc, 2) |> in_units_of(u"km"), [-180, 0],color=:black, linestyle=:solid)
    end


    ft = dominant_facies!(ax_wheel, header, data, show=:both)
    pos_min = header.axes.x |> in_units_of(u"km") |> minimum
    pos_max = header.axes.x |> in_units_of(u"km") |> maximum
    for st_bdr in st_div
        lines!(ax_wheel, [pos_min, pos_max], fill(st_bdr, 2) |> in_units_of(u"Myr"),
        color = :black,
        linestyle=:dash)
    end
    for (i, st_lab) in enumerate(st_labels)
        text!(ax_wheel, 15u"km" |> in_units_of(u"km"), st_pos[i] |> in_units_of(u"Myr"); text = st_lab,
        align=(:center, :center) )
    end
    for loc in sampl_loc
        lines!(ax_wheel, fill(loc, 2) |> in_units_of(u"km"), [header.axes.t[1], header.axes.t[end]] |> in_units_of(u"Myr"),
        color=pos_color, linewidth=pos_linewidth,
        linestyle=pos_linestyle)
    end

    Colorbar(fig[3,1], ft, vertical = false, ticks=(1:3, ["Euphotic", "Oligophotic", "Aphotic"]), label="Dominant Facies")

    min_sl = minimum(header.sea_level) |> in_units_of(u"m")
    max_sl = maximum(header.sea_level) |> in_units_of(u"m")
    fig_ext = 4
    text_adjust = 2
    ax_sl = Axis(fig[2, 2], title = "Sea level curve",xlabel = "Eustatic sea level [m]",
    limits = ((min_sl, max_sl + fig_ext), (header.axes.t[1] |> in_units_of(u"Myr"), header.axes.t[end]|> in_units_of(u"Myr"))))
    lines!(ax_sl, header.sea_level |> in_units_of(u"m") , header.axes.t|> in_units_of(u"Myr"))
    for st_bdr in st_div
        lines!(ax_sl, [min_sl, max_sl + fig_ext], fill(st_bdr, 2) |> in_units_of(u"Myr"),
        color = :black, linestyle = :dash)
    end
    for (i, st_lab) in enumerate(st_labels)
        text!(ax_sl, max_sl + text_adjust, st_pos[i] |> in_units_of(u"Myr"); text = st_lab ,
        align = (:center, :center))
    end

    linkyaxes!(ax_wheel, ax_sl)
    Label(fig[1,1, TopLeft()], "A", fontsize=subfig_label_fontsize)
    Label(fig[2,1, TopLeft()], "B", fontsize=subfig_label_fontsize)
    Label(fig[2,2, TopLeft()], "C", fontsize=subfig_label_fontsize)
    save(name, fig)
end

plot_sfig3() = plot_supp_fig( TAG1, "figs/sm/sfig3_platform_detail.png")
plot_sfig4() = plot_supp_fig(TAG2, "figs/sm/sfig4_ramp_detail.png")


plot_fig1 = function()
    fig = Figure(size=(1200, 1000), backgroundcolor=:gray80)

    header, data = read_slice("data/$(TAG1).h5", :profile)
    ax = Axis(fig[1:2, 1:2])
    sediment_profile!(ax, header, data, show_unconformities = true, show_coeval_lines = false, title = "Platform profile")
    coeval_lines!(ax, header, data, age_depth_model(header, data), [0.25:0.5:3.75;]u"Myr", linewidth = 2, color = :black, linestyle = :solid)
    ax.xlabel = distance_label
    ax.ylabel = depth_label
    Label(fig[1:2, 1:2, TopLeft()], "A", fontsize = subfig_label_fontsize)

    header, data = read_slice("data/$(TAG2).h5", :profile)
    ax = Axis(fig[1:2, 4:5])
    sediment_profile!(ax, header, data, show_unconformities = true, show_coeval_lines = false, title = "Ramp profile")
    coeval_lines!(ax, header, data, age_depth_model(header, data), [0.25:0.5:3.75;]u"Myr", linewidth = 2, color = :black, linestyle = :solid)
    ax.xlabel = distance_label
    ax.ylabel = depth_label
    Label(fig[1:2, 4:5, TopLeft()], "D",fontsize = subfig_label_fontsize)

    ax_prod1 = Axis(fig[1, 3])
    max_depth = minimum(header.initial_topography)
    h5open("data/$(TAG1).h5", "r") do fid
        production_curve!(ax_prod1, fid["input"], max_depth=max_depth)
    end
    ax_prod1.ylabel = water_depth_label
    ax_prod1.xlabel = production_label
    ax_prod1.title = "Production Platform"
    Label(fig[1, 3, TopLeft()], "B", fontsize = subfig_label_fontsize)

    ax_prod2 = Axis(fig[2, 3])
    max_depth = minimum(header.initial_topography)
    h5open("data/$(TAG2).h5", "r") do fid
        production_curve!(ax_prod2, fid["input"], max_depth=max_depth)
    end
    ax_prod2.ylabel = water_depth_label
    ax_prod2.xlabel = production_label
    ax_prod2.title = "Production Ramp"
    linkxaxes!(ax_prod1, ax_prod2)

    Label(fig[2,3, TopLeft()], "C",fontsize = subfig_label_fontsize)

    ax = Axis(fig[3, 4:5])
    header, data = read_slice("data/$(TAG2).h5", :profile)
    dominant_facies!(ax, header, data)
    Label(fig[3,4:5, TopLeft()], "G",fontsize = subfig_label_fontsize)

    ax = Axis(fig[3, 1:2])
    header, data = read_slice("data/$(TAG1).h5", :profile)
    dominant_facies!(ax, header, data)
    Label(fig[3, 1:3, TopLeft()], "E",fontsize = subfig_label_fontsize)

    ax = Axis(fig[3, 3], title = "Eustatic sea level",xlabel = "Sea level [m]",
    limits = (nothing, (header.axes.t[1] |> in_units_of(u"Myr"), header.axes.t[end]|> in_units_of(u"Myr"))))
    lines!(ax, header.sea_level |> in_units_of(u"m") , header.axes.t|> in_units_of(u"Myr"))
    ax.ylabel = time_axis_label
    Label(fig[3,3, TopLeft()], "F",fontsize = subfig_label_fontsize)

    save("figs/ms/fig1.png", fig)


end

function plot_sfig1()
    fig = summary_plot("data/$(TAG1)_prerun.h5")
    save("figs/sm/sfig1_$(TAG1)_prerun_summary.png", fig)
end

function plot_sfig2()
    fig = summary_plot("data/$(TAG2)_prerun.h5")
    save("figs/sm/sfig2_$(TAG2)_prerun_summary.png", fig)  
end

function plot_sfig20()
    fig = summary_plot("data/platform_to_ramp.h5")    
    save("figs/sm/sfig20_platform_to_ramp.png", fig) 
end

function plot_sfig21()
    fig = summary_plot("data/ramp_to_platform.h5")
    save("figs/sm/sfig21_ramp_to_platform.png", fig)
end

end

Plot.plot_fig1()

Plot.plot_sfig1()
Plot.plot_sfig2()
Plot.plot_sfig3()
Plot.plot_sfig4()
Plot.plot_sfig20()
Plot.plot_sfig21()