module Extract_data

using CarboKitten
using Unitful
using CarboKitten.Export: data_export, CSV

const TAG = "example_Anna"
const PATH = "data"

function main()
data_export(
        CSV(tuple.(10:20:70, 25),
            :age_depth_model => "$(PATH)/$(TAG)_adm.csv",
            :water_depth => "$(PATH)/$(TAG)_wd.csv",
            :metadata => "$(PATH)/$(TAG).toml"),
        "$(PATH)/$(TAG).h5")
end

end

Extract_data.main()