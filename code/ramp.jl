using Unitful
using CarboKitten
using CarboKitten.Export: read_slice, data_export, read_volume
using CarboKitten.Utility: in_units_of
using Random

using CarboKitten.Denudation: Dissolution, EmpiricalDenudation, PhysicalErosion

using CarboKitten.Models: WithDenudation as WDn

Random.seed!(255827)

const PATH = "data"
const TAG = "ramp"

const FACIES_RAMP = [
    WDn.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate=100u"m/Myr",
            extinction_coefficient=0.4u"m^-1",
            saturation_intensity=60u"W/m^2"
        ),
        transport_coefficient=40.0u"m/yr",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.0023u"m/yr",
        name = "euphotic"),
    WDn.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate=300u"m/Myr",
            extinction_coefficient=0.1u"m^-1",
            saturation_intensity=60u"W/m^2"
        ),
        transport_coefficient=40.0u"m/yr",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.0023u"m/yr",
        name = "oligophotic"),
    WDn.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        production = BenthicProduction(
            maximum_growth_rate=50u"m/Myr",
            extinction_coefficient=0.005u"m^-1",
            saturation_intensity=60u"W/m^2"
        ),
        transport_coefficient=30u"m/yr",
        reactive_surface=10u"m^2/m^3",
        mass_density=2730u"kg/m^3",
        infiltration_coefficient=0.5,
        erodibility = 0.0023u"m/yr",
        name = "aphotic")
]

include("constants.jl")

const INPUT_INIT_RUN = WDn.Input(
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
    lithification_time = LITHIFICATION_TIME_RAMP,
    denudation = DENUDATION_EMP)

println("Running ramp pre-run")
run_model(Model{WDn}, INPUT_INIT_RUN, "$(PATH)/$(TAG)_prerun.h5")

header, volume = read_volume("$(PATH)/$(TAG)_prerun.h5", :topography)
init_matrix = get_init_topography(header, volume)

const INPUT_MAIN_RUN = WDn.Input(
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
    lithification_time = LITHIFICATION_TIME_RAMP,
    denudation = DENUDATION_EMP)
    
println("Running ramp main run")
run_model(Model{WDn}, INPUT_MAIN_RUN, "$(PATH)/$(TAG).h5")
