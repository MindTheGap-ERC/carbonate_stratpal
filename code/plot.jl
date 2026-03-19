
module Plot_Profile 

using CairoMakie
using CarboKitten.Visualization: sediment_profile!, summary_plot, coeval_lines!
using CarboKitten.Export: read_slice
using Unitful

TAG1 = "platform"
TAG2 = "ramp"

function plot_platform()
    header, data = read_slice("data/$(TAG1).h5", :profile)
    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1])
    sediment_profile!(ax, header, data, show_unconformities = true, show_coeval_lines = false)
    #coeval_lines!(ax, header, data, [0.25:0.5:3.75;]u"Myr", linewidth = 3, color = :black, linestyle = :solid)
    save("figs/$(TAG1).png", fig)
    fig = summary_plot("data/$(TAG1).h5")
    save("figs/$(TAG1)_summary.png", fig)
end

function plot_ramp()
    header, data = read_slice("data/$(TAG2).h5", :profile)
    fig = Figure(size=(1000, 600))
    ax = Axis(fig[1, 1])
    sediment_profile!(ax, header, data, show_unconformities = true, show_coeval_lines = false)
    #coeval_lines!(ax, header, data, [0.25:0.5:3.75;]u"Myr", linewidth = 3, color = :black, linestyle = :solid)
    save("figs/$(TAG2).png", fig)
    fig = summary_plot("data/$(TAG2).h5")
    save("figs/$(TAG2)_summary.png", fig)
end


function plot_preruns()
    fig = summary_plot("data/$(TAG1)_prerun.h5")
    save("figs/$(TAG1)_prerun_summary.png", fig)
    fig = summary_plot("data/$(TAG2)_prerun.h5")
    save("figs/$(TAG2)_prerun_summary.png", fig)   
end
end

Plot_Profile.plot_platform()
Plot_Profile.plot_ramp()
Plot_Profile.plot_preruns()