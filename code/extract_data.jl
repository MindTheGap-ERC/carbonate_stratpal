module Extract_data

using CarboKitten
using Unitful
using CarboKitten.Export: read_slice, data_export, CSV

const PATH = "data"
const TAG1 = "platform"
const TAG2 = "ramp"

function main()
header, profile = read_slice("$(PATH)/$(TAG1).h5", :profile)
    columns = [profile[i] for i in 11:20:71]
    data_export(
        CSV(:sediment_accumulation_curve => "$(PATH)/$(TAG1)_sac.csv",
            :age_depth_model => "$(PATH)/$(TAG1)_adm.csv",
            :stratigraphic_column => "$(PATH)/$(TAG1)_sc.csv",
            :water_depth => "$(PATH)/$(TAG1)_wd.csv",
            :metadata => "$(PATH)/$(TAG1).toml"),
         header,
         columns)
header, profile = read_slice("$(PATH)/$(TAG2).h5", :profile)
    columns = [profile[i] for i in 10:20:70]
    data_export(
        CSV(:sediment_accumulation_curve => "$(PATH)/$(TAG2)_sac.csv",
            :age_depth_model => "$(PATH)/$(TAG2)_adm.csv",
            :stratigraphic_column => "$(PATH)/$(TAG2)_sc.csv",
            :water_depth => "$(PATH)/$(TAG2)_wd.csv",
            :metadata => "$(PATH)/$(TAG2).toml"),
         header,
         columns)
end
end

Extract_data.main()