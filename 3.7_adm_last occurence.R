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

#to destruct hiatus use: , lwd_destr = 0
#use approxfun for linnear interpolation

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

#Changes of system tracts
LST_to_TST = 0.5      # Lowstand system tract to Transgressive system tract
TST_to_HST = 1        # Transgressive system tract to Highstand system tract
HST_to_FSST = 1.5     # Highstand system tract to Falling Stage system tract

### Fossil abundance at adm_4km
p3(rate = 500, from = min_time(adm_4km), to = max_time(adm_4km)) |>    # constant rate in time domain
  time_to_strat(adm_4km, destructive = FALSE) |>                       # transform into depth domain
  hist(main = "Last occurrence at 4 km with constant extinction rate",
       sub = "system tracts seperated by vertical lines",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100),
       col="white")
abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="red",lwd=2,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="red",lwd=2,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="red",lwd=2,lty='dashed')
strat_to_time(7.5,adm_4km)

### Fossil abundance at adm_9km
p3(rate = 500, from = min_time(adm_9km), to = max_time(adm_9km)) |>    # constant rate in time domain
  time_to_strat(adm_9km, destructive = FALSE) |>                       # transform into depth domain
  hist(main = "Last occurrence at 9 km with constant extinction rate",
       sub = "system tracts seperated by vertical lines",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100),
       col="white")
abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="red",lwd=2,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="red",lwd=2,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="red",lwd=2,lty='dashed')

### Fossil abundance at adm_12km
p3(rate = 500, from = min_time(adm_12km), to = max_time(adm_12km)) |>  # constant rate in time domain
  time_to_strat(adm_12km, destructive = FALSE) |>                      # transform into depth domain
  hist(main = "Last occurrence at 12 km with constant extinction rate ",
       sub = "system tracts seperated by vertical lines",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100),
       col="white")
abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="red",lwd=2,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="red",lwd=2,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="red",lwd=2,lty='dashed')


### Fossil abundance at adm_15km
p3(rate = 500, from = min_time(adm_15km), to = max_time(adm_15km)) |>  # constant rate in time domain
  time_to_strat(adm_15km, destructive = FALSE) |>                      # transform into depth domain
  hist(main = "Last occurrence at 15 km with constant extinction rate",
       sub = "system tracts seperated by vertical lines",
       ylab = "Last occurrence",
       xlab = "Stratigraphic height [m]",
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100),
       col="white")
abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="red",lwd=2,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="red",lwd=2,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="red",lwd=2,lty='dashed')


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
abline(v=1.75,col="black",lwd=2,lty='dashed')
abline(v=2.25,col="black",lwd=2,lty='dashed')
abline(v=2.75,col="black",lwd=2,lty='dashed')
abline(v=(0.32+1.25),col="red",lwd=2,lty='dashed')

############### Scenario 1 Lowstand System Tract
### time: 1.25 - 1.75 Myr   #shifted: 0 - 0.5 Myr
#4km
h1.4 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
time_to_strat(0,adm_4km)                                    # start of extinction in depth domain
time_to_strat(0.25,adm_4km,destructive=FALSE)               # peak of extinction in depth domain
time_to_strat(0.5,adm_4km,destructive=FALSE)                # end of extinction in depth domain
ext1.4 <- cut(h1.4$breaks, c(-Inf,20.16365,Inf))            # separate extinction
plot(h1.4,main = "Last occurrence of taxa at 4 km",          # histogram decoration
     sub = "Extinction during the Lowstand System Tract",
     ylab = "Last occurrence",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.4])                          # extinction is red, 
  abline(v=6.228328,col="black",lwd=2,lty='dashed')      # line at peak of extinction                                                  # from peak-0.2Myr to peak+0.2Myr 
plot(adm_4km)

#9km
h1.9 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))
time_to_strat(0,adm_9km,destructive=FALSE)
time_to_strat(0.25,adm_9km)
time_to_strat(0.5,adm_9km,destructive=FALSE)
ext1.9 <- cut(h1.9$breaks, c(-Inf,60.40385,Inf))
plot(h1.9,main = "Last occurrence of taxa at 9 km",
     sub = "Extinction during the Lowstand System Tract",
     ylab = "Last occurrence",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.9])
  abline(v=5.515268,col="black",lwd=2,lty='dashed')
plot(adm_9km)

#12km
h1.12 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))
time_to_strat(0,adm_12km)
time_to_strat(0.25,adm_12km,destructive=FALSE)
time_to_strat(0.5,adm_12km)
ext1.12 <- cut(h1.12$breaks, c(-Inf,10.34304,Inf))
plot(h1.12,main = "Last occurrence of taxa at 12 km",
     sub = "Extinction during the Lowstand System Tract",
     ylab = "Last occurrence",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.12])
  abline(v=6.542699,col="black",lwd=2,lty='dashed')
plot(adm_12km)

#15km
h1.15 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_15km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))
time_to_strat(0,adm_15km)
time_to_strat(0.25,adm_15km,destructive=FALSE)
time_to_strat(0.5,adm_15km,destructive=FALSE)
ext1.15 <- cut(h1.15$breaks, c(-Inf,9.168291,Inf))
plot(h1.15,main= "Last occurrence of taxa at 15 km",
     sub = "Extinction during the Lowstand system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.15])
  abline(v=6.837498,col="black",lwd=2,lty='dashed')
plot(adm_15km)

############### Scenario 2 Transgressive System Tract
### time: 1.75 - 2.25 Myr   # Shifted: 0.5 - 1 Myr
#4km
h2.4 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
time_to_strat(0.5,adm_4km,destructive=FALSE)
time_to_strat(0.75,adm_4km)
time_to_strat(1,adm_4km)
ext2.4 <- cut(h2.4$breaks, c(-Inf,20.16365,69.7831,Inf))
plot(h2.4,main= "Last occurrence of taxa at 4 km",
     sub = "Extinction during the Transgressive system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.4])
  abline(v=45.33873,col="black",lwd=2,lty='dashed')

#9km
h2.9 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))
time_to_strat(0.5,adm_9km,destructive=FALSE)
time_to_strat(0.75,adm_9km,destructive=FALSE)
time_to_strat(1,adm_9km,destructive=FALSE)
ext2.9 <- cut(h2.9$breaks, c(-Inf,60.40385,123.6621,Inf))
plot(h2.9,main= "Last occurrence of taxa at 9 km",
     sub = "Extinction during the Transgressive system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.9])
  abline(v=89.86919,col="black",lwd=2,lty='dashed')

#12km
h2.12 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))
time_to_strat(0.5,adm_12km)
time_to_strat(0.75,adm_12km)
time_to_strat(1,adm_12km)
ext2.12 <- cut(h2.12$breaks, c(-Inf,10.34304,16.90142,Inf))
plot(h2.12,main= "Last occurrence of taxa at 12 km",
     sub = "Extinction during the Transgressive system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.12])
  abline(v=14.16978,col="black",lwd=2,lty='dashed')

#15km
h2.15 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))
time_to_strat(0.5,adm_15km,destructive=FALSE)
time_to_strat(0.75,adm_15km,destructive=FALSE)
time_to_strat(1,adm_15km,destructive=FALSE)
ext2.15 <- cut(h2.15$breaks, c(-Inf,9.168291,18.54824,Inf))
plot(h2.15,main= "Last occurrence of taxa at 15 km",
     sub = "Extinction during the Transgressive system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.15])
  abline(v=14.28246,col="black",lwd=2,lty='dashed')

############### Scenario 3 Highstand System Tract
### time: 2.25 - 2.75 Myr     # Shifted: 1 - 1.5 Myr
#4km
h3.4 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
time_to_strat(1,adm_4km)
time_to_strat(1.25,adm_4km)
time_to_strat(1.5,adm_4km,destructive=FALSE)
ext3.4 <- cut(h3.4$breaks, c(-Inf,69.7831,94.22824,Inf))
plot(h3.4,main= "Last occurrence of taxa at 4 km",
     sub = "Extinction during the Highstand system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.4])
  abline(v=86.14923,col="black",lwd=2,lty='dashed')

#9km
h3.9 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))
time_to_strat(1,adm_9km,destructive=FALSE)
time_to_strat(1.25,adm_9km,destructive=FALSE)
time_to_strat(1.5,adm_9km,destructive=FALSE)
ext3.9 <- cut(h3.9$breaks, c(-Inf,123.6621,150.853,Inf))
plot(h3.9,main= "Last occurrence of taxa at 9 km",
     sub = "Extinction during the Highstand system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.9])
abline(v=142.6458,col="black",lwd=2,lty='dashed')

#12km
h3.12 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))
time_to_strat(1,adm_12km)
time_to_strat(1.25,adm_12km)
time_to_strat(1.5,adm_12km)
ext3.12 <- cut(h3.12$breaks, c(-Inf,16.90142,25.97025,Inf))
plot(h3.12,main= "Last occurrence of taxa at 12 km",
     sub = "Extinction during the Highstand system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.12])
abline(v=20.85556,col="black",lwd=2,lty='dashed')

#15km
h3.15 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))
time_to_strat(1,adm_15km,destructive=FALSE)
time_to_strat(1.25,adm_15km)
time_to_strat(1.5,adm_15km,destructive=FALSE)
ext3.15 <- cut(h3.15$breaks, c(-Inf,18.54824,24.37772,Inf))
plot(h3.15,main= "Last occurrence of taxa at 15 km",
     sub = "Extinction during the Highstand system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.15])
abline(v=19.69094,col="black",lwd=2,lty='dashed')

############### Scenario 4 Falling Stage System Tract
### time: 2.75 - 3.25 Myr       # Shifted: 1.5 - 2 Myr
#4km
h4.4 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
time_to_strat(1.5,adm_4km,destructive=FALSE)
time_to_strat(1.75,adm_4km,destructive=FALSE)
time_to_strat(2,adm_4km,destructive=FALSE)
ext4.4 <- cut(h4.4$breaks, c(-Inf,94.22824,100.0478,Inf))
plot(h4.4,main= "Last occurrence of taxa at 4 km",
     sub = "Extinction during the Falling Stage system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.4])
  abline(v=97.68081,col="black",lwd=2,lty='dashed')

#9km
h4.9 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))
time_to_strat(1.5,adm_9km,destructive=FALSE)
time_to_strat(1.75,adm_9km,destructive=FALSE)
time_to_strat(2,adm_9km)
ext4.9 <- cut(h4.9$breaks, c(-Inf,150.853,156.665,Inf))
plot(h4.9,main= "Last occurrence of taxa at 9 km",
     sub = "Extinction during the Falling Stage system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.9])
  abline(v=154.2088,col="black",lwd=2,lty='dashed')

#12km
h4.12 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))
time_to_strat(1.5,adm_12km,destructive=FALSE)
time_to_strat(1.75,adm_12km)
time_to_strat(2,adm_12km)
ext4.12 <- cut(h4.12$breaks, c(-Inf,25.97025,94.43012,Inf))
plot(h4.12,main= "Last occurrence of taxa at 12 km",
     sub = "Extinction during the Falling Stage system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.12])
  abline(v=33.81947,col="black",lwd=2,lty='dashed')

#15km
h4.15 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))
time_to_strat(1.5,adm_15km,destructive=FALSE)
time_to_strat(1.75,adm_15km,destructive=FALSE)
time_to_strat(2,adm_15km,destructive=FALSE)
ext4.15 <- cut(h4.15$breaks, c(-Inf,24.37772,34.35685,Inf))
plot(h4.15,main= "Last occurrence of taxa at 15 km",
     sub = "Extinction during the Falling Stage system tract",
     ylab = "Last occurrence frequency",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.15])
  abline(v=29.49508,col="black",lwd=2,lty='dashed')

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
