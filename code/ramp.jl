module Ramp

using Unitful
using CarboKitten

const PATH = "data"

const TAG = "ramp"

const FACIES = [
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
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=40.0u"m/yr"),
    ALCAP.Facies(
        viability_range = (4, 10),
        activation_range = (6, 10),
        maximum_growth_rate=50u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=30u"m/yr")
]

const PERIOD1 = 2.0u"Myr"
const AMPLITUDE1 = 8.0u"m"
const PERIOD2 = 0.2u"Myr"
const AMPLITUDE2 = 2.0u"m"

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(100, 50), phys_scale=150.0u"m"),
    time=TimeProperties(
        Δt=100u"yr",
        steps=10000),
    output=Dict(
        :profile => OutputSpec(slice=(:, 25), write_interval=1)),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> AMPLITUDE1 * sin(2π * t / PERIOD1) + AMPLITUDE2 * sin(2π * t / PERIOD2),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=50,
    depositional_resolution=500.0u"m",
    facies=FACIES)


function main()
    run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")

end

end

Ramp.main()