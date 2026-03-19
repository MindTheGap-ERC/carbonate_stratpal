const PERIOD1 = 2.0u"Myr"
const AMPLITUDE1 = 15.0u"m"
const PERIOD2 = 0.2u"Myr"
const AMPLITUDE2 = 2.5u"m"
const DELTA_T = 100u"yr"
const STEPS_INIT = 10000
const STEPS_MAIN = 40000
const PHYS_SCALE = 150u"m"
const GRID_SIZE_Y = 20
const GRID_SIZE_X = 150
const BOX = Box{Coast}(grid_size=(GRID_SIZE_X, GRID_SIZE_Y), phys_scale=PHYS_SCALE)
const SUBSIDENCE_RATE = 30u"m/Myr"
const DISINTGRATION_RATE = 50.0u"m/Myr"
const INSOLATION = 400.0u"W/m^2"
const SEDIMENT_BUFFER_SIZE = 50
const DEPOSITIONAL_RESOLUTION = 500.0u"m"
const LITHIFICATION_TIME_PLATFORM = 100u"yr"
const LITHIFICATION_TIME_RAMP = 1u"kyr"
const INIT_SL = t -> 0u"m"
const MAIN_SL = t -> AMPLITUDE1 * sin(2π * t / PERIOD1) + AMPLITUDE2 * sin(2π * t / PERIOD2)

const DENUDATION_EMP = EmpiricalDenudation(precip = 1.0u"m")

function get_init_topography(header, volume)    
    n_steps = size(volume.sediment_thickness)[3]
    bedrock = header.initial_topography .- header.axes.t[end] * header.subsidence_rate
    grid_size = (length(header.axes.x), length(header.axes.y))
    result = Array{Float64, 2}(undef, grid_size...)
    result[:, :] = (volume.sediment_thickness[:,:,n_steps] .+ bedrock) |> in_units_of(u"m")
end