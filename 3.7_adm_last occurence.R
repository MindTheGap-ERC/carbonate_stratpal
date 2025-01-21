data_kitten = read.csv("C:/Users/sidne/OneDrive/Documenten/AA Utrecht/Guided research/Git repositories/Stratigraphic Paleobiology for Phenotypic Evolution/CarboKitten.jl/data/output/alcap-example3_adm.csv")
#grid: 20km x 50km
#time: 4 Myr
#Sea level from Holland & Patzkowsky 2015

library(admtools)
library(StratPal)
# PC file name: "3.7"
#ADM's: 500m:500m:15km
t = data_kitten$time..Myr.
t_mod = t[ t >= 1.25 & t <= 3.25]                          # Subset data
t_mod_shifted = t_mod - min(t_mod)
###         ###  4 km
h4 = data_kitten$adm8..m.                                  # 4km
h4_mod = h4[ t >= 1.25 & t <= 3.25]
h4_mod_shifted = h4_mod - min(h4_mod)

adm4_mod = tp_to_adm(t = t_mod, h4_mod)                    # adm modified to show 1.25-3.25 Myr
adm_4km = tp_to_adm(t = t_mod_shifted, h4_mod_shifted)     # adm modified & shifted
plot(adm4_mod)
plot(adm_4km)
###         ###  9 km
h9 = data_kitten$adm18..m.                                  # 9km
h9_mod = h9[ t >= 1.25 & t <= 3.25]
h9_mod_shifted = h9_mod - min(h9_mod)

adm9_mod = tp_to_adm(t = t_mod, h9_mod)                    # adm modified to show 1.25-3.25 Myr
adm_9km = tp_to_adm(t = t_mod_shifted, h9_mod_shifted)     # adm modified & shifted
plot(adm_9km)
###         ###  12 km
h12 = data_kitten$adm24..m.                                # 12km
h12_mod = h12[ t >= 1.25 & t <= 3.25]
h12_mod_shifted = h12_mod - min(h12_mod)

adm12_mod = tp_to_adm(t = t_mod, h12_mod)                  # adm modified to show 1.25-3.25 Myr
adm_12km = tp_to_adm(t = t_mod_shifted, h12_mod_shifted)   # adm modified & shifted
plot(adm_12km)
###         ###  15 km
h15 = data_kitten$adm30..m.                                # 15km
h15_mod = h15[ t >= 1.25 & t <= 3.25]
h15_mod_shifted = h15_mod - min(h15_mod)

adm15_mod = tp_to_adm(t = t_mod, h15_mod)                  # adm modified to show 1.25-3.25 Myr
adm_15km = tp_to_adm(t = t_mod_shifted, h15_mod_shifted)   # adm modified & shifted
plot(adm_15km)
### example of previous adm code:
## adm_6km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm12..m.)
###


# plotting
plot(adm_1km); text(2,190,"age depth model at 1km")
plot(adm_2km); text(2,190,"age depth model at 2km")
plot(adm_3km); text(2,190,"age depth model at 3km")
plot(adm_4km); text(2,190,"age depth model at 4km")
hist(get_hiat_duration(adm_4km))
#longer gaps should correspond to 4th order SL change
plot(adm_5km); text(2,200,"age depth model at 5km")
hist(get_hiat_duration(adm_5km))
plot(adm_6km); text(2,200,"age depth model at 6km")
hist(get_hiat_duration(adm_6km))
plot(adm_7km); text(2,200,"age depth model at 7km")
hist(get_hiat_duration(adm_7km))
plot(adm_8km); text(2,210,"age depth model at 8km")
hist(get_hiat_duration(adm_8km))
plot(adm_9km); text(2,210,"age depth model at 9km")
hist(get_hiat_duration(adm_9km))
plot(adm_10km); text(2,210,"age depth model at 10km")
hist(get_hiat_duration(adm_10km))
plot(adm_11km); text(2,210,"age depth model at 11km")
hist(get_hiat_duration(adm_11km))
plot(adm_12km); text(2,210,"age depth model at 12km")
hist(get_hiat_duration(adm_12km))
plot(adm_13km); text(2,200,"age depth model at 13km")
hist(get_hiat_duration(adm_13km))
plot(adm_14km); text(2,90,"age depth model at 14km")
hist(get_hiat_duration(adm_14km))
plot(adm_15km); text(2,65,"age depth model at 15km")
hist(get_hiat_duration(adm_15km))
plot(adm_16km); text(2,60,"age depth model at 16km")
plot(adm_17km); text(2,60,"age depth model at 17km")
plot(adm_18km); text(2,65,"age depth model at 18km")
plot(adm_19km); text(2,65,"age depth model at 19km")
plot(adm_20km); text(2,65,"age depth model at 20km")
#to destruct hiatus use: , lwd_destr = 0
#use approxfun for linnear interpolation

adm <- list(adm_1km,adm_2km,adm_3km,adm_4km,adm_5km,adm_6km,adm_7km,adm_8km,
            adm_9km,adm_10km,adm_11km,adm_12km,adm_13km,adm_14km,adm_15km,
            adm_16km,adm_17km,adm_18km,adm_19km,adm_20km)


###Plotting adm_4km off-shore
plot(adm_4km, text(0.5,90,"age depth model at 4km"))
plot_sed_rate_t(adm_4km)

###Random walk evolution at adm_4km
seq(from = min_time(adm_4km), to = max_time(adm_4km), by = 0.01) |> # sample every 10 kyr over the interval covered by the adm 
  random_walk(sigma = 1, mu = 3) |>                                 # simulate random walk
  time_to_strat(adm_4km, destructive = FALSE) |>                    # transform data from time to strat domain
  plot(type = "l",                                                  # plot
       orientation = "lr",
       xlab = paste0("Stratigraphic height [", get_L_unit(adm_4km), "]"),
       ylab = "Trait value",
       main = "Trait evolution 4 km from shore")

### Fossil abundance at adm_4km
p3(rate = 500, from = min_time(adm_4km), to = max_time(adm_4km)) |> # constant rate in time domain
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(xlab = "Stratigraphic height [m]",                           # plot
       main = "Fossil abundance 4 km offshore",
       ylab = "# Fossils",
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))

########## Eustatic Sea Level
AMP1 = 15
PER1 = 2
AMP2 = 2.5
PER2 = 0.2

t = seq(1.25, 3.25, by=0.001)
sl = (AMP1*sin(2*pi*t/PER1)) + (AMP2*sin(2*pi*t/PER2))
plot(t,sl,type='l',
     main="Eustatic Sea Level",
     xlab="time",
     ylab="meter")

############### Scenario 1 Lowstand System Tract
### time: 1.25 - 1.75 Myr   #shifted: 0 - 0.5 Myr
#4km
p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 4 km offshore",
       sub = "Extinction during the Lowstand System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
time_to_strat(0,adm_4km)
time_to_strat(0.25,adm_4km,destructive=FALSE)
time_to_strat(0.5,adm_4km,destructive=FALSE)
plot(adm_4km)

#9km
p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 9 km offshore",
       sub = "Extinction during the Lowstand System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))
time_to_strat(0,adm_9km,destructive=FALSE)
time_to_strat(0.25,adm_9km)
time_to_strat(0.5,adm_9km,destructive=FALSE)
plot(adm_9km)

#12km
p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 12 km offshore",
       sub = "Extinction during the Lowstand System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))
time_to_strat(0,adm_12km)
time_to_strat(0.25,adm_12km,destructive=FALSE)
time_to_strat(0.5,adm_12km)
plot(adm_12km)

#15km
p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_15km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 15 km offshore",
       sub = "Extinction during the Lowstand System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))
time_to_strat(0,adm_15km)
time_to_strat(0.25,adm_15km,destructive=FALSE)
time_to_strat(0.5,adm_15km,destructive=FALSE)
plot(adm_15km)

############### Scenario 2 Transgressive System Tract
### time: 1.75 - 2.25 Myr   # Shifted: 0.5 - 1 Myr
#4km
p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 4 km offshore",
       sub = "Extinction during the Transgressive System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
time_to_strat(0.5,adm_4km,destructive=FALSE)
time_to_strat(0.75,adm_4km)
time_to_strat(1,adm_4km)
plot(adm_4km)

#9km
p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 9 km offshore",
       sub = "Extinction during the Transgressive System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))
#12km
p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 12 km offshore",
       sub = "Extinction during the Transgressive System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))
#15km
p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 15 km offshore",
       sub = "Extinction during the Transgressive System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

############### Scenario 3 Highstand System Tract
### time: 2.25 - 2.75 Myr     # Shifted: 1 - 1.5 Myr
#4km
p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 4 km offshore",
       sub = "Extinction during the Highstand System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
#9km
p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 9 km offshore",
       sub = "Extinction during the Highstand System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))
#12km
p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 12 km offshore",
       sub = "Extinction during the Highstand System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))
#15km
p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 15 km offshore",
       sub = "Extinction during the Highstand System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

############### Scenario 4 Falling Stage System Tract
### time: 2.75 - 3.25 Myr       # Shifted: 1.5 - 2 Myr
#4km
p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 4 km offshore",
       sub = "Extinction during the Falling Stage System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
#9km
p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 9 km offshore",
       sub = "Extinction during the Falling Stage System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))
#12km
p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 12 km offshore",
       sub = "Extinction during the Falling Stage System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))
#15km
p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                     # transform into depth domain
  hist(main = "Fossil abundance 15 km offshore",
       sub = "Extinction during the Falling Stage System Tract",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))
#subset
t = data_kitten$time..Myr.
h = data_kitten$adm30..m.
t_mod = t[ t >= 1.25 & t <= 3.25]
h_mod = h[ t >= 1.25 & t <= 3.25]
t_mod_shifted = t_mod - min(t_mod)
h_mod_shifted = h_mod - min(h_mod)




adm_1km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm2..m.)
plot(adm_1km); text(2,190,"age depth model at 1km")

strat_to_time(25,adm_12km)
