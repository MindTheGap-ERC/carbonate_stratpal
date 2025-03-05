module Extract_data

using CarboKitten
using Unitful
using CarboKitten.Export: data_export, CSV

const TAG = "example_Anna"
const PATH = "data"

function plot_model()
data_export(
        CSV(tuple.(10:20:70, 25),
            :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
            :metadata => "$(PATH)/$(TAG).toml"),
        "$(PATH)/$(TAG).h5")
end

end

Extract_data.plot_model()