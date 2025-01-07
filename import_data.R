data_kitten = read.csv("C:/Users/sidne/OneDrive/Documenten/AA Utrecht/Guided research/Git repositories/Stratigraphic Paleobiology for Phenotypic Evolution/CarboKitten.jl/data/output/alcap-example2_adm.csv")
#grid: 20km x 50km
#time: 4 Myr
#Sea level from Holland & Patzkowsky 2015

library(admtools)
library(StratPal)
#ADM's: 500m:500m:15km
# Computer name "3.6"
adm_1km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm2..m.)
adm_2km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm4..m.)
adm_3km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm6..m.)
adm_4km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm8..m.)
adm_5km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm10..m.)
adm_6km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm12..m.)
adm_7km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm14..m.)
adm_8km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm16..m.)
adm_9km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm18..m.)
adm_10km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm20..m.)
adm_11km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm22..m.)
adm_12km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm24..m.)
adm_13km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm26..m.)
adm_14km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm28..m.)
adm_15km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm30..m.)
adm_16km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm32..m.)
adm_17km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm34..m.)
adm_18km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm36..m.)
adm_19km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm38..m.)
adm_20km = tp_to_adm(t = data_kitten$time..Myr., data_kitten$adm40..m.)

# plotting
plot(adm_1km)
plot(adm_2km)
plot(adm_3km)
plot(adm_4km,
     main = "4km",
     xlab = "time",
     ylab = "depth")

plot(adm_5km)
plot(adm_6km)
plot(adm_7km)
plot(adm_8km)
plot(adm_9km)
plot(adm_10km)
plot(adm_11km)
plot(adm_12km)
plot(adm_13km)
plot(adm_14km)
plot(adm_15km)
plot(adm_16km)
plot(adm_17km)
plot(adm_18km)
plot(adm_19km)
plot(adm_20km)
#to destruct hiatus use: , lwd_destr = 0

adm <- list(adm_1km,adm_2km,adm_3km,adm_4km,adm_5km,adm_6km,adm_7km,adm_8km,
            adm_9km,adm_10km,adm_11km,adm_12km,adm_13km,adm_14km,adm_15km,
            adm_16km,adm_17km,adm_18km,adm_19km,adm_20km)

#Incompleteness percentage per adm #################
w_incompleteness <- c()
for(incompleteness in adm){
  v_incompleteness <- get_incompleteness(incompleteness)
  w_incompleteness <- c(w_incompleteness,v_incompleteness)
}
plot(w_incompleteness,
     main = "incompleteness per adm 1 to 20km",
     xlab = "km",
     ylab = "incompleteness percentage", 
     pch = c(16))

#Completeness percentage per adm #################
w_completeness <- c()
for(completeness in adm){
  v_completeness <- get_completeness(completeness)
  w_completeness <- c(w_completeness,v_completeness)
}
plot(w_completeness,
     main = "completeness per adm 1 to 20km",
     xlab = "km",
     ylab = "percentage",
     pch = c(16))

#Number of gaps per adm ##############################
w_hiat_no <- c()
for(hiat_no in adm){
  v_hiat_no <- get_hiat_no(hiat_no)
  w_hiat_no <- c(w_hiat_no,v_hiat_no)
}
plot(w_hiat_no,
     main = "number of gaps per adm 1 to 20km",
     xlab = "km",
     ylab = "number of gaps",
     pch = c(16))

#Average duration of hiatus per adm ##################
w_hiat_dur <- c()
for(hiat_dur in adm){
  v_hiat_dur <- mean(get_hiat_duration(hiat_dur))
  w_hiat_dur <- c(w_hiat_dur,v_hiat_dur)
}
plot(w_hiat_dur,
     main = "average duration of hiatus per adm 1 to 20km",
     xlab = "km",
     ylab = "Myr",
     pch = c(16))

#Total hiatus duration per adm #######################
w_hiat_sum <- c()
for(hiat_sum in adm){
  v_hiat_sum <- sum(get_hiat_duration(hiat_sum))
  w_hiat_sum <- c(w_hiat_sum,v_hiat_sum)
}
plot(w_hiat_sum,
     main = "total duration of gaps per adm 1 to 20km",
     xlab = "km",
     ylab = "Myr",
     pch = c(16))

#Height per adm ######################################
w_height <- c()
for(height in adm){
  v_height <- get_height(height)
  w_height <- c(w_height,v_height)
}
plot(w_height,
     main = "height",
     xlab = "?",
     ylab = "meter",
     pch = c(16))    
#plot(get_height(adm_15km))

###Plotting adm_4km off-shore
plot(adm_4km); text(2,210,"age depth model at 4km")
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
p3(rate = 200, from = min_time(adm_4km), to = max_time(adm_4km)) |> # constant rate in time domain
  time_to_strat(adm_4km, destructive = TRUE) |>                     # transform into depth domain
  hist(xlab = "Stratigraphic height [m]",                           # plot
       main = "Fossil abundance 4 km offshore",
       ylab = "# Fossils",
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 20))

# Fossil abundance at other adm's
# 8km
plot(adm_8km); text(2,210,"age depth model at 8km")
plot_sed_rate_t(adm_8km)
p3(rate = 200, from = min_time(adm_8km), to = max_time(adm_8km)) |> # constant rate in time domain
  time_to_strat(adm_8km, destructive = TRUE) |>                     # transform into depth domain
  hist(xlab = "Stratigraphic height [m]",                           # plot
       main = "Fossil abundance 8 km offshore",
       ylab = "# Fossils",
       breaks = seq(from = min_height(adm_8km), to = max_height(adm_8km), length.out = 20))

# 12km
plot(adm_12km); text(2,210,"age depth model at 12km")
plot_sed_rate_t(adm_12km)
p3(rate = 200, from = min_time(adm_12km), to = max_time(adm_12km)) |> # constant rate in time domain
  time_to_strat(adm_12km, destructive = TRUE) |>                     # transform into depth domain
  hist(xlab = "Stratigraphic height [m]",                           # plot
       main = "Fossil abundance 12 km offshore",
       ylab = "# Fossils",
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 20))

# 15km
plot(adm_15km); text(2,210,"age depth model at 15km")
plot_sed_rate_t(adm_15km)
p3(rate = 200, from = min_time(adm_15km), to = max_time(adm_15km)) |> # constant rate in time domain
  time_to_strat(adm_15km, destructive = TRUE) |>                     # transform into depth domain
  hist(xlab = "Stratigraphic height [m]",                           # plot
       main = "Fossil abundance 15 km offshore",
       ylab = "# Fossils",
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 20))
