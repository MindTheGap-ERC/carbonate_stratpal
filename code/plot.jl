
module Plot_Profile 

using CairoMakie
using CarboKitten.Visualization: sediment_profile
using CarboKitten.Export: read_slice

function plot(path, name)
    header, data = read_slice(path * name, :profile)
    fig = sediment_profile(header, data)
    save("figs/$(name).png", fig)
end

end
