using Unitful
using CarboKitten
using CarboKitten.Export: read_slice, data_export, read_volume
using CarboKitten.Utility: in_units_of

const PATH = "data"
const TAG = "ramp"

const FACIES_RAMP = [
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=200u"m/Myr",
        extinction_coefficient=0.6u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=80.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=300u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=100.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=50u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=30u"m/yr")
]

include("constants.jl")

const INPUT_INIT_RUN = ALCAP.Input(
    tag="$TAG",
    box=BOX,
    time=TimeProperties(
        Δt=DELTA_T,
        steps=STEPS_INIT),
    output=Dict(
        :profile => OutputSpec(slice=(:, div(GRID_SIZE_Y, 2)), write_interval=10),
        :topography => OutputSpec(slice = (:, :), write_interval = 1000)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=INIT_SL,
    subsidence_rate=SUBSIDENCE_RATE,
    disintegration_rate=DISINTGRATION_RATE,
    insolation=INSOLATION,
    sediment_buffer_size=SEDIMENT_BUFFER_SIZE,
    depositional_resolution=DEPOSITIONAL_RESOLUTION,
    facies=FACIES_RAMP,
    cementation_time = CEMENTATION_TIME)

println("Running initial model")
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
        :profile => OutputSpec(slice=(:, div(GRID_SIZE_Y, 2)), write_interval=10),
        :topography => OutputSpec(slice = (:, :), write_interval = 100)),
    ca_interval=1,
    initial_topography=init_matrix * u"m",
    sea_level=MAIN_SL,
    subsidence_rate=SUBSIDENCE_RATE,
    disintegration_rate=DISINTGRATION_RATE,
    insolation=INSOLATION,
    sediment_buffer_size=SEDIMENT_BUFFER_SIZE,
    depositional_resolution=DEPOSITIONAL_RESOLUTION,
    facies=FACIES_RAMP,
    cementation_time = CEMENTATION_TIME)
    
println("Running main model")
run_model(Model{ALCAP}, INPUT_MAIN_RUN, "$(PATH)/$(TAG).h5")
