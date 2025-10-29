using Unitful
using HDF5
using DelimitedFiles
import CSV as OfficialCSV
import CarboKitten
using DataFrames
using CarboKitten
using CarboKitten.Export: data_export, CSV as CSVCarbo

const x_grid_size = 100
const y_grid_size = 50
const PHYS_SCALE = 150.0u"m"
const SLOPE = 300.0

PATH = "data"
TAG = "init_carbonate_stratpal_1"

const FACIES = [
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=12.5u"m/yr")
]

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=100u"yr",
        steps=10000),
    output=Dict(
        :profile => OutputSpec(slice=(:, 25), write_interval=10),
        :topography => OutputSpec(slice = (:, :), write_interval = 100)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> 0u"m",
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=500.0u"m",
    facies=FACIES)

PATH = "data/"
TAG = "init_carbonate_stratpal_1"

function extract_topography(PATH,TAG)
    h5open("$(PATH)/$(TAG).h5", "r") do fid
        disintegration = read(fid["topography/disintegration"])[1,:,:,end]
        production = read(fid["topography/production"])[1,:,:,end]
        deposition = read(fid["topography/deposition"])[1,:,:,end]
        sediment_height = read(fid["topography/sediment_thickness"])[:,:,end]
        data_dis = DataFrame(
            disintegration, :auto
        )
        data_pro = DataFrame(
            production, :auto
        )   
        data_dep = DataFrame(
           deposition, :auto
        )
        data_sed = DataFrame(
           sediment_height, :auto
        )
        return data_dis.*1.0u"m", data_pro.*1.0u"m", data_dep.*1.0u"m", data_sed.*1.0u"m"
end
end

data_dis, data_pro, data_dep, data_sed = extract_topography(PATH,TAG)

function starting_bathy()
    init = ones(x_grid_size, y_grid_size) .*1.0u"m"
    for i in CartesianIndices(init)
        init[i]   = -(i[1]-1) .* PHYS_SCALE ./ SLOPE
    end
    return init
end

function calculate_bathymetry(data,INPUT)
    Bathy = zeros(x_grid_size, y_grid_size) .*1.0u"m"
    Bathy .= starting_bathy() .+ data .- INPUT.subsidence_rate .* INPUT.time.Δt .* INPUT.time.steps
    OfficialCSV.write("$(PATH)/$(TAG).csv", DataFrame(Bathy,:auto))
end

calculate_bathymetry(data_sed,INPUT)