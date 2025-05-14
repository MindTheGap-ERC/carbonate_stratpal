module Script

using Unitful
using CarboKitten
using CarboKitten.Export: data_export, CSV

const PATH = "data"

const TAG = "lineage-example"

const FACIES = [
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=500u"m/Myr",
        extinction_coefficient=0.8u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=50u"m/yr"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=400u"m/Myr",
        extinction_coefficient=0.1u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=25u"m/yr"),
    ALCAP.Facies(
        viability_range=(4, 10),
        activation_range=(6, 10),
        maximum_growth_rate=100u"m/Myr",
        extinction_coefficient=0.005u"m^-1",
        saturation_intensity=60u"W/m^2",
        diffusion_coefficient=35u"m/yr")
]


const PERIOD1 = 2u"Myr"
const AMPLITUDE1 = 15u"m"
const PERIOD2 = 0.2u"Myr"
const AMPLITUDE2 = 2.5u"m"

const INPUT = ALCAP.Input(
    tag="$TAG",
    box=Box{Coast}(grid_size=(200, 50), phys_scale=100.0u"m"),
    time=TimeProperties(
        Δt=0.0002u"Myr",
        steps=20000,
        write_interval=100),
    ca_interval=1,
    initial_topography=(x, y) -> -x / 300.0,
    sea_level=t -> (AMPLITUDE1*sin(2*pi*t/PERIOD1)) + (AMPLITUDE2*sin(2*pi*t/PERIOD2)),
    subsidence_rate=50.0u"m/Myr",
    disintegration_rate=50.0u"m/Myr",
    insolation=400.0u"W/m^2",
    sediment_buffer_size=100,
    depositional_resolution=0.5u"m",
    facies=FACIES)
# ~/~ end

function main()
    run_model(Model{ALCAP}, INPUT, "$(PATH)/$(TAG).h5")

    data_export(
        CSV(tuple.(5:5:200, 25),
            :sediment_accumulation_curve => "$(PATH)/$(TAG)_sac.csv",
            :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
            :stratigraphic_column => "$(PATH)/$(TAG)_sc.csv",
            #:water_depth => "$(PATH)/$(TAG)_wd.csv",
            :metadata => "$(PATH)/$(TAG).toml"),
        "$(PATH)/$(TAG).h5")
end

end

Script.main()
# ~/~ end
