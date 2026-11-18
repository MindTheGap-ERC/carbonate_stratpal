module Plot
using CarboKitten.Algorithms: skeleton
import CarboKitten.Visualization: sediment_profile!, production_curve!
using CarboKitten.Export: Header, Data, DataSlice, read_header, read_volume, read_slice, group_datasets
using CarboKitten.Utility: in_units_of
using CarboKitten.Output.Abstract: stratigraphic_column
using HDF5
using Unitful
using CarboKitten.BoundaryTrait
using CarboKitten.Stencil: convolution
using CairoMakie

const na = [CartesianIndex()]

#### Print settings ####
# Same style as the figures of the phylostrat project.
# Figures are designed at final print size for a B5 dissertation (176 x 250 mm page).
# Figure sizes are given in pt (1 pt = 1/72 in), so font sizes and line widths below
# correspond to their printed size. Raster output is written at `DPI`.
const MM = 72 / 25.4                 # pt per mm
const TEXT_WIDTH = 130 * MM          # usable text width on a B5 page [pt]
const DPI = 600
const PX_PER_UNIT = DPI / 72

const LW_MAIN = 1.0                  # data lines
const LW_SEA_LEVEL = 0.6             # sea level curve (thin: high-frequency signal)
const LW_COEVAL = 0.8                # time lines / coeval lines
const LW_SAMPLING = 0.8              # sampling positions
# unconformities were drawn with linewidth = 1 in the original 1200 px wide figure.
# Scale to keep their original size relative to the figure.
const LW_UNCONFORMITY = TEXT_WIDTH / 1200
const LW_SEALEVEL_PROFILE = 2.5      # end sea level in sediment profiles
const PANEL_LABEL_SIZE = 11
const HEADER_SIZE = 9                # figure/column headers
const ST_LABEL_SIZE = 6              # systems tract labels

const PRINT_THEME = Theme(
    fontsize = 8,
    figure_padding = (2, 6, 2, 2),
    backgroundcolor = :white,
    Axis = (
        titlesize = 8.5,
        titlefont = :bold,
        titlealign = :left,
        xlabelsize = 8,
        ylabelsize = 8,
        xticklabelsize = 7,
        yticklabelsize = 7,
        spinewidth = 0.6,
        xtickwidth = 0.6,
        ytickwidth = 0.6,
        xticksize = 3,
        yticksize = 3,
        xgridvisible = false,
        ygridvisible = false,
        xlabelpadding = 2,
        ylabelpadding = 2,
    ),
    Colorbar = (
        labelsize = 8,
        ticklabelsize = 7,
        size = 8,
        spinewidth = 0.6,
        tickwidth = 0.6,
        ticksize = 3,
    ),
    Lines = (linewidth = LW_MAIN,),
)

panel_label!(pos, txt) = Label(pos, txt; fontsize = PANEL_LABEL_SIZE, font = :bold,
    padding = (0, 6, 4, 0), halign = :right)

header_label!(pos, txt) = Label(pos, txt; fontsize = HEADER_SIZE, font = :bold, tellwidth = false)

save_print(fname, fig) = save(fname, fig; px_per_unit = PX_PER_UNIT)

#### Constants ####
const sampl_loc = [3, 7.5, 10.5, 12, 18]u"km"
const st_div = 0.25u"Myr":0.5u"Myr":3.75u"Myr"
const st_labels = ["TST1", "HST1", "RST1", "LST1", "TST2", "HST2", "RST2", "LST2", "TST3"]
const st_pos = [0.125u"Myr", 0.5u"Myr":0.5u"Myr":3.5u"Myr"..., 3.875u"Myr"]
const st_label_x = 15u"km"           # horizontal position of systems tract labels in Wheeler diagrams
const Rate = typeof(1.0u"m/Myr")
const Amount = typeof(1.0u"m")
const Length = typeof(1.0u"m")
const Time = typeof(1.0u"Myr")

const time_axis_label = "Elapsed model time [Myr]"
const water_depth_label = "Water depth [m]"
const production_label = "Production [m/Myr]"
const distance_label = "Distance from shore [km]"
const depth_label = "Depth [m]"
const sea_level_label = "Sea level [m]"

const FACIES_TICKS = (1:3, ["Euphotic", "Oligophotic", "Aphotic"])
const FACIES_LABEL = "Dominant carbonate factory"

TAG1 = "platform"
TAG2 = "ramp"

#### Low level plotting ####
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
    hiatus = skeleton(combined_acc .< 0.000001, minwidth=minwidth)

    if !isempty(hiatus[1])
        verts = [(x[pt[1]], h[pt...] |> in_units_of(u"m")) for pt in hiatus[1]]
        linesegments!(ax, vec(permutedims(verts[hiatus[2]])); kwargs...)
    end
end

function profile_plot!(ax::Axis, header::Header, data::DataSlice; color::AbstractArray, mesh_args...)
    x = header.axes.x |> in_units_of(u"km")

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

    bedrock = (initial_topography .- total_subsidence) |> in_units_of(u"m")
    lower_limit = minimum(bedrock) - 20
    band!(ax, x, lower_limit, bedrock; color=:gray, label="initial topography")
    lines!(ax, x, bedrock; color=:black, label="initial topography")
    ylims!(ax, lower_limit + 10, nothing)
    xlims!(ax, x[1], x[end])
    ax.xlabel = distance_label
    ax.ylabel = depth_label

    c = reshape(color, n_x * n_t)
    mesh!(ax, v, f; color=vcat(c, c), mesh_args...)
end

function plot_sealevel!(ax::Axis, header::Header)
    sea_level = header.sea_level[end] |> in_units_of(u"m")
    hlines!(ax, sea_level, color=:lightblue, linewidth=LW_SEALEVEL_PROFILE, label="end sea level")
end

function sediment_profile!(ax::Axis, header::Header, data::DataSlice;
                           show_unconformities::Union{Nothing,Bool,Int} = true,
                           show_coeval_lines::Union{Bool,Tuple{Int, Int},Vector{Int},Vector{Time}} = true,
                           show_sealevel::Bool = true,
                           title = "Sediment profile")
    n_facies = size(data.production)[1]

    if show_sealevel
        plot_sealevel!(ax, header)
    end

    plot = profile_plot!(argmax, ax, header, data; alpha=1.0,
        colormap=cgrad(Makie.wong_colors()[1:n_facies], n_facies, categorical=true))

    minwidth = 10
    if show_unconformities
        plot_unconformities(ax, header, data, minwidth; label = "unconformities",
                        color=:white, linestyle=:dash, linewidth=LW_UNCONFORMITY)
    end

    if show_coeval_lines
        coeval_lines!(ax, header, data, age_depth_model(header, data), time_lines_in_range(header, st_div);
            color=:black, linestyle=:dash, linewidth=LW_COEVAL)
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

    ax.ylabel = time_axis_label
    ax.xlabel = distance_label
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
    n_facies, n_x, n_t = size(data.production)
    total_subsidence = (header.axes.t[end] - header.axes.t[1]) * header.subsidence_rate
    initial_topography = header.initial_topography[data.slice...]
    sc = stratigraphic_column(data)
    h = repeat(initial_topography .- total_subsidence, 1, n_t+1)
    @views h[:, 2:end] .+= cumsum(sum(sc, dims=1)[1,:,:], dims=2)
    return h
end

# time lines that fall within the simulated time interval (e.g. for the 1 Myr pre-runs)
time_lines_in_range(header::Header, time_lines) =
    Time[t for t in time_lines if header.axes.t[1] < t < header.axes.t[end]]

#### Panels ####
function my_profile!(ax::Axis, header::Header, data::DataSlice;
                     time_lines = st_div,
                     sampling_positions = nothing,
                     title = "Sediment profile")
    sediment_profile!(ax, header, data; show_coeval_lines = false, show_unconformities = true)
    coeval_lines!(ax, header, data, age_depth_model(header, data), time_lines_in_range(header, time_lines);
        linewidth = LW_COEVAL, color = :black, linestyle = :dash)
    if !isnothing(sampling_positions)
        vlines!(ax, sampling_positions |> in_units_of(u"km");
            linewidth = LW_SAMPLING, linestyle = :solid, color = :black)
    end
    ax.xlabel = distance_label
    ax.ylabel = depth_label
    ax.title = title
    return ax
end

function my_wheeler!(ax::Axis, header::Header, data::DataSlice;
                     time_lines = st_div,
                     sampling_positions = nothing,
                     show_st_labels = false,
                     title = "Wheeler diagram")
    ft = dominant_facies!(ax, header, data, show = :both)
    tl = time_lines_in_range(header, time_lines)
    if !isempty(tl)
        hlines!(ax, tl |> in_units_of(u"Myr"); linewidth = LW_COEVAL, color = :black, linestyle = :dash)
    end
    if !isnothing(sampling_positions)
        vlines!(ax, sampling_positions |> in_units_of(u"km");
            linewidth = LW_SAMPLING, linestyle = :solid, color = :black)
    end
    if show_st_labels
        for (i, st_lab) in enumerate(st_labels)
            text!(ax, st_label_x |> in_units_of(u"km"), st_pos[i] |> in_units_of(u"Myr");
                text = st_lab, align = (:center, :center), fontsize = ST_LABEL_SIZE)
        end
    end
    ax.xlabel = distance_label
    ax.ylabel = time_axis_label
    ax.title = title
    return ft
end

function sea_level_panel!(ax::Axis, header::Header;
                          time_lines = st_div,
                          show_st_labels = false)
    sl = header.sea_level |> in_units_of(u"m")
    t = header.axes.t |> in_units_of(u"Myr")
    lines!(ax, sl, t; linewidth = LW_SEA_LEVEL, color = :blue)
    tl = time_lines_in_range(header, time_lines)
    if !isempty(tl)
        hlines!(ax, tl |> in_units_of(u"Myr"); linewidth = LW_COEVAL, color = :black, linestyle = :dash)
    end
    min_sl, max_sl = extrema(sl)
    if min_sl == max_sl # constant sea level (pre-runs)
        min_sl, max_sl = min_sl - 1, max_sl + 1
    end
    if show_st_labels
        # extend the axis to the right to make room for the labels
        text_x = max_sl + 0.15 * (max_sl - min_sl)
        for (i, st_lab) in enumerate(st_labels)
            text!(ax, text_x, st_pos[i] |> in_units_of(u"Myr");
                text = st_lab, align = (:center, :center), fontsize = ST_LABEL_SIZE)
        end
        max_sl = max_sl + 0.3 * (max_sl - min_sl)
    end
    limits!(ax, min_sl, max_sl, t[1], t[end])
    ax.title = "Eustatic sea level"
    ax.xlabel = sea_level_label
    ax.ylabel = time_axis_label
    ax.xticks = LinearTicks(3)
    return ax
end

function production_panel!(ax::Axis, input::String, max_depth; title = "Production profile")
    h5open(input, "r") do fid
        production_curve!(ax, fid["input"], max_depth = max_depth)
    end
    ax.xlabel = production_label
    ax.ylabel = water_depth_label
    ax.title = title # set after plotting, CarboKitten overwrites it
    return ax
end

facies_colorbar!(pos, ft) = Colorbar(pos, ft; vertical = false, flipaxis = false,
    ticks = FACIES_TICKS, label = FACIES_LABEL,
    width = Relative(0.6), tellheight = true)

#### Figures ####
function plot_fig1()
    fig = Figure(size = (TEXT_WIDTH, 165 * MM))

    header_pl, data_pl = read_slice("data/$(TAG1).h5", :profile)
    header_ra, data_ra = read_slice("data/$(TAG2).h5", :profile)

    ax_profile_pl = Axis(fig[1:2, 1:2])
    my_profile!(ax_profile_pl, header_pl, data_pl, title = "Platform profile")
    panel_label!(fig[1:2, 1:2, TopLeft()], "A")

    ax_profile_ra = Axis(fig[1:2, 4:5])
    my_profile!(ax_profile_ra, header_ra, data_ra, title = "Ramp profile")
    panel_label!(fig[1:2, 4:5, TopLeft()], "D")

    max_depth = minimum(header_ra.initial_topography)
    ax_prod_pl = Axis(fig[1, 3])
    production_panel!(ax_prod_pl, "data/$(TAG1).h5", max_depth, title = "Platform\nproduction")
    panel_label!(fig[1, 3, TopLeft()], "B")

    ax_prod_ra = Axis(fig[2, 3])
    production_panel!(ax_prod_ra, "data/$(TAG2).h5", max_depth, title = "Ramp\nproduction")
    panel_label!(fig[2, 3, TopLeft()], "C")
    linkxaxes!(ax_prod_pl, ax_prod_ra)
    # narrow column: fewer ticks to avoid overlapping tick labels
    ax_prod_pl.xticks = LinearTicks(3)
    ax_prod_ra.xticks = LinearTicks(3)

    ax_wheeler_pl = Axis(fig[3, 1:2])
    ft = my_wheeler!(ax_wheeler_pl, header_pl, data_pl, title = "Platform Wheeler diagram")
    panel_label!(fig[3, 1:2, TopLeft()], "E")

    ax_sl = Axis(fig[3, 3])
    sea_level_panel!(ax_sl, header_pl)
    panel_label!(fig[3, 3, TopLeft()], "F")

    ax_wheeler_ra = Axis(fig[3, 4:5])
    my_wheeler!(ax_wheeler_ra, header_ra, data_ra, title = "Ramp Wheeler diagram")
    panel_label!(fig[3, 4:5, TopLeft()], "G")

    linkxaxes!(ax_profile_pl, ax_wheeler_pl)
    linkxaxes!(ax_profile_ra, ax_wheeler_ra)
    linkyaxes!(ax_wheeler_pl, ax_sl, ax_wheeler_ra)
    # shared time axis: drop redundant labels on the sea level panel
    hideydecorations!(ax_sl, ticks = false, grid = false)

    facies_colorbar!(fig[4, 1:5], ft)

    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 8)
    colgap!(fig.layout, 3, 16) # room for the production x label next to the ramp profile y label
    save_print("figs/ms/fig1.png", fig)
end

function plot_supp_fig(tag, name, case)
    header, data = read_slice("data/$(tag).h5", :profile)

    fig = Figure(size = (TEXT_WIDTH, 125 * MM))
    header_label!(fig[0, 1:2], case)

    ax_profile = Axis(fig[1, 1])
    my_profile!(ax_profile, header, data, sampling_positions = sampl_loc)
    panel_label!(fig[1, 1, TopLeft()], "A")

    ax_wheeler = Axis(fig[2, 1])
    ft = my_wheeler!(ax_wheeler, header, data, sampling_positions = sampl_loc, show_st_labels = true)
    panel_label!(fig[2, 1, TopLeft()], "B")

    ax_sl = Axis(fig[2, 2])
    sea_level_panel!(ax_sl, header, show_st_labels = true)
    panel_label!(fig[2, 2, TopLeft()], "C")

    linkxaxes!(ax_profile, ax_wheeler)
    linkyaxes!(ax_wheeler, ax_sl)
    hideydecorations!(ax_sl, ticks = false, grid = false)

    facies_colorbar!(fig[3, 1], ft)

    colsize!(fig.layout, 1, Relative(0.72))
    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 8)
    save_print(name, fig)
end

plot_sfig3() = plot_supp_fig(TAG1, "figs/sm/sfig3_platform_detail.png", "Platform")
plot_sfig4() = plot_supp_fig(TAG2, "figs/sm/sfig4_ramp_detail.png", "Ramp")

# replaces CarboKitten's `summary_plot`, same layout as the summary plots of phylostrat
function my_summary_plot(input::String, fig_name::String, case::String; time_lines = st_div)
    header, data = read_slice(input, :profile)

    fig = Figure(size = (TEXT_WIDTH, 125 * MM))
    header_label!(fig[0, 1:2], case)

    ax_profile = Axis(fig[1, 1])
    my_profile!(ax_profile, header, data, time_lines = time_lines)
    panel_label!(fig[1, 1, TopLeft()], "A")

    ax_production = Axis(fig[1, 2])
    production_panel!(ax_production, input, minimum(header.initial_topography))
    panel_label!(fig[1, 2, TopLeft()], "B")

    ax_wheeler = Axis(fig[2, 1])
    ft = my_wheeler!(ax_wheeler, header, data, time_lines = time_lines)
    panel_label!(fig[2, 1, TopLeft()], "C")

    ax_sl = Axis(fig[2, 2])
    sea_level_panel!(ax_sl, header, time_lines = time_lines)
    panel_label!(fig[2, 2, TopLeft()], "D")

    linkxaxes!(ax_profile, ax_wheeler)
    linkyaxes!(ax_wheeler, ax_sl)
    hideydecorations!(ax_sl, ticks = false, grid = false)

    facies_colorbar!(fig[3, 1], ft)

    colsize!(fig.layout, 1, Relative(0.68))
    rowgap!(fig.layout, 6)
    colgap!(fig.layout, 8)
    save_print(fig_name, fig)
end

# pre-runs have constant sea level, so no systems tracts: use evenly spaced time lines
const prerun_time_lines = 0.25u"Myr":0.25u"Myr":0.75u"Myr"
plot_sfig1() = my_summary_plot("data/$(TAG1)_prerun.h5", "figs/sm/sfig1_$(TAG1)_prerun_summary.png", "Platform pre-run";
    time_lines = prerun_time_lines)
plot_sfig2() = my_summary_plot("data/$(TAG2)_prerun.h5", "figs/sm/sfig2_$(TAG2)_prerun_summary.png", "Ramp pre-run";
    time_lines = prerun_time_lines)
plot_sfig20() = my_summary_plot("data/platform_to_ramp.h5", "figs/sm/sfig20_platform_to_ramp.png", "Platform to ramp")
plot_sfig21() = my_summary_plot("data/ramp_to_platform.h5", "figs/sm/sfig21_ramp_to_platform.png", "Ramp to platform")

set_theme!(PRINT_THEME)

end

Plot.plot_fig1()

Plot.plot_sfig1()
Plot.plot_sfig2()
Plot.plot_sfig3()
Plot.plot_sfig4()
Plot.plot_sfig20()
Plot.plot_sfig21()
