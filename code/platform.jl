
using Unitful
using CarboKitten
using CarboKitten.Export: read_slice, data_export, read_volume
using CarboKitten.Utility: in_units_of


const TAG = "platform"
const PATH = "data"
const FACIES_PLATFORM = [
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

const PERIOD1 = 2.0u"Myr"
const AMPLITUDE1 = 15.0u"m"
const PERIOD2 = 0.2u"Myr"
const AMPLITUDE2 = 2.5u"m"
const DELTA_T = 100u"yr"
const STEPS_INIT = 10000
const STEPS_MAIN = 40000
const PHYS_SCALE = 150u"m"
const BOX = Box{Coast}(grid_size=(150, 50), phys_scale=PHYS_SCALE)
const SUBSIDENCE_RATE = 30u"m/Myr"

function get_init_topography(header, volume)    
    n_steps = size(volume.sediment_thickness)[3]
    bedrock = header.initial_topography .- header.axes.t[end] * header.subsidence_rate
    grid_size = (length(header.axes.x), length(header.axes.y))
    result = Array{Float64, 2}(undef, grid_size...)
    result[:, :] = (volume.sediment_thickness[:,:,n_steps] .+ bedrock) |> in_units_of(u"m")
end






const INPUT_INIT_RUN = ALCAP.Input(
    tag="$TAG",
    box=BOX,
    time=TimeProperties(
        Δt=DELTA_T,
        steps=STEPS_INIT),
    output=Dict(
        :profile => OutputSpec(slice=(:, 25), write_interval=10),
        :topography => OutputSpec(slice = (:, :), write_interval = 1000)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> 0u"m",
    subsidence_rate=SUBSIDENCE_RATE,
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=500.0u"m",
    facies=FACIES_PLATFORM,
    cementation_time = 1u"kyr")

println("Running initial model\n")
run_model(Model{ALCAP}, INPUT_INIT_RUN, "$(PATH)/$(TAG)_prerun.h5")

header, volume = read_volume("$(PATH)/$(TAG)_prerun.h5", :topography)

init_matrix = get_init_topography(header, volume)

const INPUT_MAIN_RUN = ALCAP.Input(
    tag="$TAG",
    box=BOX,
    time=TimeProperties(
        Δt=DELTA_T,
        steps=STEPS_MAIN),
    output=Dict(
        :profile => OutputSpec(slice=(:, 25), write_interval=10),
        :topography => OutputSpec(slice = (:, :), write_interval = 100)),
    ca_interval=1,
    initial_topography=init_matrix * u"m", #initial_topography1(path),#(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE1 * sin(2π * t / PERIOD1) + AMPLITUDE2 * sin(2π * t / PERIOD2),
    subsidence_rate=SUBSIDENCE_RATE,
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=500.0u"m",
    facies=FACIES_PLATFORM,
    cementation_time = 1u"kyr")

println("Running main model\n")
run_model(Model{ALCAP}, INPUT_MAIN_RUN, "$(PATH)/$(TAG).h5")

