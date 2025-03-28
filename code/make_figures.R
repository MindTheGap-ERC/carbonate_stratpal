data_kitten = read.csv("C:/Users/sidne/OneDrive/Documenten/AA Utrecht/Guided research/Git repositories/Stratigraphic Paleobiology for Phenotypic Evolution/lineage_ecol/data/alcap-example7_adm.csv")
water_depth = read.csv("C:/Users/sidne/OneDrive/Documenten/AA Utrecht/Guided research/Git repositories/Stratigraphic Paleobiology for Phenotypic Evolution/lineage_ecol/data/alcap-example7_wd.csv")

#grid: 20km x 50km
#time: 4 Myr
#Sea level from Holland & Patzkowsky 2015

library(admtools)
library(StratPal)

### Eustatic Sea Level curve
AMP1 = 15                           # Amplitude and period of 3rd order fluctuation
PER1 = 2                            
AMP2 = 2.5                          # Amplitude and period of 4th order fluctuation
PER2 = 0.2
t = seq(1.25, 3.25, by=0.001)       # time frame
sl = (AMP1*sin(2*pi*t/PER1)) + (AMP2*sin(2*pi*t/PER2))    # the curve
plot(t_mod_shifted,sl,type='l',     # Plot curve
     main="Eustatic Sea Level",
     xlab="time [Myr]",
     ylab="[m]")
abline(v=0.5,col="cyan4",lwd=3,lty='dashed')    # separate different system tracts
abline(v=1,col="cyan4",lwd=3,lty='dashed')
abline(v=1.5,col="cyan4",lwd=3,lty='dashed')
text(x=0.2, y=0, "Lowstand 
     system tract",cex=0.9,col='darkblue',font=3)
text(x=0.75, y=12, "Transgressive 
     system tract",cex=0.9,col='darkblue',font=3)
text(x=1.2, y=0, "Highstand 
     system tract",cex=0.9,col='darkblue',font=3)
text(x=1.75, y=12, "Falling stage 
     system tract",cex=0.9,col='darkblue',font=3)

# 40 age depth models. From 500m to 20km, every 500m
t = data_kitten$time..Myr.
t_mod = t[t >= 1.25 & t <= 3.25]                           # Subset data
t_mod_shifted = t_mod - min(t_mod)
#Changes of system tracts
LST_to_TST = 0.5      # Lowstand system tract to Transgressive system tract
TST_to_HST = 1        # Transgressive system tract to Highstand system tract
HST_to_FSST = 1.5     # Highstand system tract to Falling Stage system tract

###  4 km     adm
h4 = data_kitten$adm8..m.                           
h4_mod = h4[t >= 1.25 & t <= 3.25]                         # modify to show 1.25-3.25 Myr
h4_mod_shifted = h4_mod - min(h4_mod)                       # shift to 0-2 Myr
adm4_mod = tp_to_adm(t = t_mod, h4_mod)                    
adm_4km = tp_to_adm(t = t_mod_shifted, h4_mod_shifted)     
plot(h4,type='l')
plot(adm_4km,lwd_acc = 2,lwd_destr = 0)                     # plot                               
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty="dashed")
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty="dashed")
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty="dashed")
title(main = "age depth model 4km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
text(x=0.21, y=40, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.75, y=80, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.2, y=20, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.8, y=60, "FSST",cex=0.8,col='darkblue',font=3)

w4 = water_depth$wd8..m.                                   # 4km water depth (wd)        
w4_mod = w4[t >= 1.25 & t <= 3.25]                         # modify and shift
w4_mod_shifted = w4_mod - min(w4_mod)
wd4_mod = tp_to_adm(t = t_mod, w4_mod)                    
wd_4km = tp_to_adm(t = t_mod_shifted, w4_mod_shifted)       
plot(w4_mod)                                                # plot
plot(w4) 
plot(w4,type='l')
title(main="water depth at 4km",
      ylab="[m]",
      xlab="time [Myr]")
plot(w4,type='l',main="complete wd 4km",
     xlab="time [Myr]",ylab="[m]")

###  9 km       adm
h9 = data_kitten$adm18..m.                                 
h9_mod = h9[ t >= 1.25 & t <= 3.25]                         # modify to show 1.25-3.25 Myr
h9_mod_shifted = h9_mod - min(h9_mod)                       # shift to 0-2 Myr
adm9_mod = tp_to_adm(t = t_mod, h9_mod)                     
adm_9km = tp_to_adm(t = t_mod_shifted, h9_mod_shifted)      
plot(h9,type='l')
plot(adm_9km,lwd_acc = 2,lwd_destr = 0)                     # plot
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty="dashed")
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty="dashed")
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty="dashed")
title(main = "age depth model 9km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
text(x=0.22, y=50, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.78, y=50, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.28, y=50, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.8, y=50, "FSST",cex=0.8,col='darkblue',font=3)

w9 = water_depth$wd18..m.                                   # 9km water depth (wd)                                
w9_mod = w9[ t >= 1.25 & t <= 3.25]                         # modify and shift
w9_mod_shifted = w9_mod - min(w9_mod)
wd9_mod = tp_to_adm(t = t_mod, w9_mod)                    
wd_9km = tp_to_adm(t = t_mod_shifted, w9_mod_shifted)
plot(wd_9km)                                                # plot
title(main="water depth at 9km",
      ylab="[m]",
      xlab="time [Myr]")
plot(w9,type='l',main="complete wd 9km",
     xlab="time [Myr]",ylab="[m]")


### 10km adm
h10 = data_kitten$adm20..m.                                 # 10km adm to compare to 12km
h10_mod = h10[ t >= 1.25 & t <= 3.25]
h10_mod_shifted = h10_mod - min(h10_mod)
adm10_mod = tp_to_adm(t = t_mod, h10_mod)                  
adm_10km = tp_to_adm(t = t_mod_shifted, h10_mod_shifted)   
plot(adm_10km,lwd_acc = 2,lwd_destr = 0)
title('10km adm')
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')
### 11km adm
h11 = data_kitten$adm22..m.                                 # 11km adm to compare to 12km
h11_mod = h11[ t >= 1.25 & t <= 3.25]
h11_mod_shifted = h11_mod - min(h11_mod)
adm11_mod = tp_to_adm(t = t_mod, h11_mod)                  
adm_11km = tp_to_adm(t = t_mod_shifted, h11_mod_shifted)   
plot(adm_11km,lwd_acc = 2,lwd_destr = 0)
title('11km adm')
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')

###  12 km    adm
h12 = data_kitten$adm24..m.                                
h12_mod = h12[ t >= 1.25 & t <= 3.25]                       # modify to show 1.25-3.25 Myr
h12_mod_shifted = h12_mod - min(h12_mod)                    # shift to 0-2 Myr
adm12_mod = tp_to_adm(t = t_mod, h12_mod)           
adm_12km = tp_to_adm(t = t_mod_shifted, h12_mod_shifted) 

plot(adm_12km,lwd_acc = 2,lwd_destr = 0)                    # plot
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')
title(main = "age depth model 12km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
text(x=0.21, y=40, "Lowstand 
     system tract",cex=0.6,col='darkblue',font=3)
text(x=0.7, y=50, "Transgressive 
     system tract",cex=0.6,col='darkblue',font=3)
text(x=1.2, y=60, "Highstand
     system tract",cex=0.6,col='darkblue',font=3)
text(x=1.7, y=70, "Falling stage 
     system tract",cex=0.6,col='darkblue',font=3)

w12 = water_depth$wd24..m.                                  # 12km water depth (wd)                                
w12_mod = w12[ t >= 1.25 & t <= 3.25]                       # modify and shift
w12_mod_shifted = w12_mod - min(w12_mod)
wd12_mod = tp_to_adm(t = t_mod, w12_mod)                    
wd_12km = tp_to_adm(t = t_mod_shifted, w12_mod_shifted)     
plot(wd_12km)                                               # plot
title(main="water depth at 12km",
      ylab="[m]",
      xlab="time [Myr]")
plot(w12,type='l',main="complete wd 12km",
     xlab="time [Myr]",ylab="[m]")

### 13 km adm
h13 = data_kitten$adm26..m.                                 # 13 km adm to compare to 12km                           
h13_mod = h13[ t >= 1.25 & t <= 3.25]
h13_mod_shifted = h13_mod - min(h13_mod)
adm13_mod = tp_to_adm(t = t_mod, h13_mod)                   
adm_13km = tp_to_adm(t = t_mod_shifted, h13_mod_shifted)    
plot(adm_13km,lwd_acc = 2,lwd_destr = 0)                                              
title('13km adm')
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')

### 15 km    adm
h15 = data_kitten$adm30..m.                                 # 15km adm
h15_mod = h15[ t >= 1.25 & t <= 3.25]                       # modified to show 1.25-3.25 Myr
h15_mod_shifted = h15_mod - min(h15_mod)                    # shift to 0-2 Myr
adm15_mod = tp_to_adm(t = t_mod, h15_mod)                  
adm_15km = tp_to_adm(t = t_mod_shifted, h15_mod_shifted)   

plot(adm_15km,lwd_acc = 2,lwd_destr = 0)                    # plot
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')
title(main = "age depth model 15km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
text(x=0.2, y=15, "Lowstand 
     system tract",cex=0.8,col='darkblue',font=3)
text(x=0.72, y=5, "Transgressive 
     system tract",cex=0.8,col='darkblue',font=3)
text(x=1.2, y=10, "Highstand
     system tract",cex=0.8,col='darkblue',font=3)
text(x=1.7, y=15, "Falling stage 
     system tract",cex=0.8,col='darkblue',font=3)

w15 = water_depth$wd30..m.                                    # 15km water depth (wd)                               
w15_mod = w15[ t >= 1.25 & t <= 3.25]                         # modify and shift
w15_mod_shifted = w15_mod - min(w15_mod)
wd15_mod = tp_to_adm(t = t_mod, w15_mod)                    
wd_15km = tp_to_adm(t = t_mod_shifted, w15_mod_shifted)     
par(new=TRUE)
plot(wd_15km)                                                 # plot
title(main="water depth at 15 km",
      ylab="[m]",
      xlab="time [Myr]")
plot(w15,type='l',main="complete wd 15km",
     xlab="time [Myr]")

###############################################################################
### 4km Last occurrences        constant extinction rate
p3(rate = 500, from = min_time(adm_4km), to = max_time(adm_4km)) |>    # constant rate in time domain
  time_to_strat(adm_4km, destructive = FALSE) |>                       # transform into depth domain
  hist(main = "Last occurrence at 4 km with constant extinction rate",
       ylab = "Last occurrence per meter",
       xlab = "Stratigraphic height [m]",
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100),
       col="white")
abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=8, y=95, "LST",cex=0.9,col='darkblue',font=3)
text(x=45, y=95, "TST",cex=0.9,col='darkblue',font=3)
text(x=82, y=95, "HST",cex=0.9,col='darkblue',font=3)
text(x=95, y=95, "FSST",cex=0.9,col='darkblue',font=3)

### 9km Last ccurrences         constant extinction rate
p3(rate = 500, from = min_time(adm_9km), to = max_time(adm_9km)) |>    # constant rate in time domain
  time_to_strat(adm_9km, destructive = FALSE) |>                       # transform into depth domain
  hist(main = "Last occurrence at 9 km with constant extinction rate",
       sub = "system tracts seperated by vertical lines",
       ylab = "Last occurrence per meter",
       xlab = "Stratigraphic height [m]",
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100),
       col="white")
abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=25, y=95, "LST",cex=0.9,col='darkblue',font=3)
text(x=90, y=95, "TST",cex=0.9,col='darkblue',font=3)
text(x=135, y=95, "HST",cex=0.9,col='darkblue',font=3)
text(x=157, y=95, "FSST",cex=0.9,col='darkblue',font=3)

### 12km Last occurrence        constant extinction rate
p3(rate = 500, from = min_time(adm_12km), to = max_time(adm_12km)) |>  # constant rate in time domain
  time_to_strat(adm_12km, destructive = FALSE) |>                      # transform into depth domain
  hist(main = "Last occurrence at 12 km with constant extinction rate ",
       sub = "system tracts seperated by vertical lines",
       ylab = "Last occurrence per meter",
       xlab = "Stratigraphic height [m]",
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100),
       col="white")
abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=79, "LST",cex=0.9,col='darkblue',font=3)
text(x=13, y=79, "TST",cex=0.9,col='darkblue',font=3)
text(x=21, y=79, "HST",cex=0.9,col='darkblue',font=3)
text(x=40, y=79, "FSST",cex=0.9,col='darkblue',font=3)
(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE))
max_height(adm_12km)

### 15km Last occurrence          constant extinction rate
p3(rate = 500, from = min_time(adm_15km), to = max_time(adm_15km)) |>  # constant rate in time domain
  time_to_strat(adm_15km, destructive = FALSE) |>                      # transform into depth domain
  hist(main = "Last occurrence at 15 km with constant extinction rate",
       sub = "system tracts seperated by vertical lines",
       ylab = "Last occurrence per meter",
       xlab = "Stratigraphic height [m]",
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100),
       col="white")
abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=100, "LST",cex=0.9,col='darkblue',font=3)
text(x=13, y=100, "TST",cex=0.9,col='darkblue',font=3)
text(x=22, y=100, "HST",cex=0.9,col='darkblue',font=3)
text(x=30, y=100, "FSST",cex=0.9,col='darkblue',font=3)
par(new=TRUE)
plot(wd_15km)                                                 # plot
axis(4)
mtext("[m]",side=4)
title(main="water depth at 15 km",
      ylab="[m]",
      xlab="time [Myr]")
plot(w15,type='l',main="complete wd 15km",
     xlab="time [Myr]",ylab="[m]")
axis(4)
mtext("[m]",side=4)
################################################################################
### Extinction event simulations
# Shape of extinction
tp = c(0,0.5,0.75,1,2)                            # points in time
ext_rate = c(1,1,25,1,1)                          # rate of extinction
plot(tp,ext_rate,type='l',lwd=3,                  # graph simply showing the extinction
     main="extinction event",
     xlab="time [Myr]",
     ylab="extinction rate [species/m]")

### Scenario 1 - Lowstand System Tract extincion - time: 1.25-1.75 Myr (shifted to 0-0.5 Myr)
#4km LST
h1.4 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))

ext1.4 <- cut(h1.4$breaks, c(-Inf,(time_to_strat(0.5,adm_4km,destructive=FALSE)),Inf))  # distinguishing the timing of extinction event
plot(h1.4,main = "Last occurrence of taxa at 4 km",                  # histogram with coloring
     sub = "Mass extinction during Lowstand system tract ",
     font.sub=3,
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.4])                          # extinction is red, rest is white 
  abline(v=(time_to_strat(0.25,adm_4km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=10, y=145, "LST",cex=0.9,col='darkblue',font=3)
text(x=45, y=145, "TST",cex=0.9,col='darkblue',font=3)
text(x=82, y=145, "HST",cex=0.9,col='darkblue',font=3)
text(x=100, y=145, "FSST",cex=0.9,col='darkblue',font=3)

#9km LST
h1.9 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))

ext1.9 <- cut(h1.9$breaks, c(-Inf,(time_to_strat(0.5,adm_9km,destructive=FALSE)),Inf))  # distinguishing the timing of extinction event
plot(h1.9,main = "Last occurrence of taxa at 9 km",                  # histogram with coloring
     sub = "Mass extinction during the Lowstand System Tract",
     font.sub=3,
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.9])
  abline(v=(time_to_strat(0.25,adm_9km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=25, y=100, "LST",cex=0.9,col='darkblue',font=3)
text(x=90, y=100, "TST",cex=0.9,col='darkblue',font=3)
text(x=138, y=100, "HST",cex=0.9,col='darkblue',font=3)
text(x=158, y=100, "FSST",cex=0.9,col='darkblue',font=3)

#12km LST
h1.12 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences  
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))

ext1.12 <- cut(h1.12$breaks, c(-Inf,(time_to_strat(0.5,adm_12km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h1.12,main = "Last occurrence of taxa at 12 km",                # histogram with coloring
     sub = "Mass extinction during the Lowstand System Tract",
     font.sub=3,
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.12])
  abline(v=(time_to_strat(0.25,adm_12km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')  
  abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=3.8, y=70, "LST",cex=0.9,col='darkblue',font=3)
text(x=13.5, y=70, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=70, "HST",cex=0.9,col='darkblue',font=3)
text(x=50, y=70, "FSST",cex=0.9,col='darkblue',font=3)

#15km LST
h1.15 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |> 
  time_to_strat(adm_15km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences  
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

ext1.15 <- cut(h1.15$breaks, c(-Inf,(time_to_strat(0.5,adm_15km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h1.15,main= "Last occurrence of taxa at 15 km",                 # histogram with coloring
     sub = "Mass extinction during the Lowstand system tract",
     font.sub=3,
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("red","white")[ext1.15])
  abline(v=(time_to_strat(0.25,adm_15km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=110, "LST",cex=0.9,col='darkblue',font=3)
text(x=14, y=110, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=110, "HST",cex=0.9,col='darkblue',font=3)
text(x=30, y=110, "FSST",cex=0.9,col='darkblue',font=3)

### Scenario 2 - Transgressive System Tract extinction - time: 1.75-2.25 Myr (shifted to 0.5-1 Myr)
#4km TST
h2.4 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))

ext2.4 <- cut(h2.4$breaks, c(-Inf,(time_to_strat(0.5,adm_4km,destructive=FALSE)),(time_to_strat(1,adm_4km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h2.4,main= "Last occurrence of taxa at 4 km",                   # histogram with coloring
     sub = "Mass extinction during the Transgressive system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.4])
  abline(v=(time_to_strat(0.75,adm_4km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=10, y=45, "LST",cex=0.9,col='darkblue',font=3)
text(x=45, y=45, "TST",cex=0.9,col='darkblue',font=3)
text(x=82, y=45, "HST",cex=0.9,col='darkblue',font=3)
text(x=100, y=45, "FSST",cex=0.9,col='darkblue',font=3)

#9km TST
h2.9 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))

ext2.9 <- cut(h2.9$breaks, c(-Inf,(time_to_strat(0.5,adm_9km,destructive=FALSE)),(time_to_strat(1,adm_9km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h2.9,main= "Last occurrence of taxa at 9 km",                   # histogram with coloring
     sub = "Mass extinction during the Transgressive system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.9])
  abline(v=(time_to_strat(0.75,adm_9km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate 
  abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=25, y=90, "LST",cex=0.9,col='darkblue',font=3)
text(x=100, y=90, "TST",cex=0.9,col='darkblue',font=3)
text(x=138, y=90, "HST",cex=0.9,col='darkblue',font=3)
text(x=157, y=90, "FSST",cex=0.9,col='darkblue',font=3)

#12km TST
h2.12 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))

ext2.12 <- cut(h2.12$breaks, c(-Inf,(time_to_strat(0.5,adm_12km,destructive=FALSE)),(time_to_strat(1,adm_12km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h2.12,main= "Last occurrence of taxa at 12 km",                 # histogram with coloring
     sub = "Mass extinction during the Transgressive system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.12])
  abline(v=(time_to_strat(0.75,adm_12km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=190, "LST",cex=0.9,col='darkblue',font=3)
text(x=13.5, y=190, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=190, "HST",cex=0.9,col='darkblue',font=3)
text(x=50, y=190, "FSST",cex=0.9,col='darkblue',font=3)

#15km TST
h2.15 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

ext2.15 <- cut(h2.15$breaks, c(-Inf,(time_to_strat(0.5,adm_15km,destructive=FALSE)),(time_to_strat(1,adm_15km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h2.15,main= "Last occurrence of taxa at 15 km",                 # histogram with coloring
     sub = "Mass extinction during the Transgressive system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext2.15])
  abline(v=(time_to_strat(0.75,adm_15km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=85, "LST",cex=0.9,col='darkblue',font=3)
text(x=14, y=85, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=85, "HST",cex=0.9,col='darkblue',font=3)
text(x=30, y=85, "FSST",cex=0.9,col='darkblue',font=3)

#### Scenario 3 - Highstand System Tract extinction - time: 2.25-2.75 Myr (shifted to 1-1.5 Myr)
#4km HST
h3.4 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))
time_to_strat(1,adm_4km)
time_to_strat(1.25,adm_4km)
time_to_strat(1.5,adm_4km,destructive=FALSE)
ext3.4 <- cut(h3.4$breaks, c(-Inf,(time_to_strat(1,adm_4km,destructive=FALSE)),(time_to_strat(1.5,adm_4km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h3.4,main= "Last occurrence of taxa at 4 km",                   # histogram with coloring
     sub = "Mass extinction during the Highstand system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.4])
  abline(v=(time_to_strat(1.25,adm_4km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=10, y=145, "LST",cex=0.9,col='darkblue',font=3)
text(x=45, y=145, "TST",cex=0.9,col='darkblue',font=3)
text(x=81, y=145, "HST",cex=0.9,col='darkblue',font=3)
text(x=100, y=145, "FSST",cex=0.9,col='darkblue',font=3)

#9km HST
h3.9 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))

ext3.9 <- cut(h3.9$breaks, c(-Inf,(time_to_strat(1,adm_9km,destructive=FALSE)),(time_to_strat(1.5,adm_9km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h3.9,main= "Last occurrence of taxa at 9 km",                  # histogram with coloring
     sub = "Mass extinction during the Highstand system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.9])
  abline(v=(time_to_strat(1.25,adm_9km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=30, y=153, "LST",cex=0.9,col='darkblue',font=3)
text(x=92, y=153, "TST",cex=0.9,col='darkblue',font=3)
text(x=140, y=153, "HST",cex=0.9,col='darkblue',font=3)
text(x=157, y=153, "FSST",cex=0.9,col='darkblue',font=3)

#12km HST
h3.12 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))

ext3.12 <- cut(h3.12$breaks, c(-Inf,(time_to_strat(1,adm_12km,destructive=FALSE)),(time_to_strat(1.5,adm_12km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h3.12,main= "Last occurrence of taxa at 12 km",                 # histogram with coloring
     sub = "Mass extinction during the Highstand system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.12])
  abline(v=(time_to_strat(1.25,adm_12km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=61, "LST",cex=0.9,col='darkblue',font=3)
text(x=14, y=61, "TST",cex=0.9,col='darkblue',font=3)
text(x=22, y=61, "HST",cex=0.9,col='darkblue',font=3)
text(x=55, y=61, "FSST",cex=0.9,col='darkblue',font=3)

#15km HST
h3.15 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

ext3.15 <- cut(h3.15$breaks, c(-Inf,(time_to_strat(1,adm_15km,destructive=FALSE)),(time_to_strat(1.5,adm_15km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h3.15,main= "Last occurrence of taxa at 15 km",                 # histogram with coloring
     sub = "Mass extinction during the Highstand system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext3.15])
  abline(v=(time_to_strat(1.25,adm_15km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=115, "LST",cex=0.9,col='darkblue',font=3)
text(x=14, y=115, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=115, "HST",cex=0.9,col='darkblue',font=3)
text(x=30, y=115, "FSST",cex=0.9,col='darkblue',font=3)

### Scenario 4 - Falling Stage System Tract extinction - time: 2.75-3.25 Myr (shifted to 1.5-2 Myr)
#4km FSST
h4.4 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_4km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_4km), to = max_height(adm_4km), length.out = 100))

ext4.4 <- cut(h4.4$breaks, c(-Inf,(time_to_strat(1.5,adm_4km,destructive=FALSE)),(time_to_strat(2,adm_4km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h4.4,main= "Last occurrence of taxa at 4 km",                   # histogram with coloring
     sub = "Mass extinction during the Falling Stage system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.4])
  abline(v=(time_to_strat(1.75,adm_4km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=10, y=230, "LST",cex=0.9,col='darkblue',font=3)
text(x=45, y=230, "TST",cex=0.9,col='darkblue',font=3)
text(x=82, y=230, "HST",cex=0.9,col='darkblue',font=3)
text(x=98, y=230, "FSST",cex=0.9,col='darkblue',font=3)

#9km FSST
h4.9 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_9km, destructive = FALSE) |>                     # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_9km), to = max_height(adm_9km), length.out = 100))

ext4.9 <- cut(h4.9$breaks, c(-Inf,(time_to_strat(1.5,adm_9km,destructive=FALSE)),(time_to_strat(2,adm_9km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h4.9,main= "Last occurrence of taxa at 9 km",                   # histogram with coloring
     sub = "Mass extinction during the Falling Stage system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.9])
  abline(v=(time_to_strat(1.75,adm_9km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=25, y=240, "LST",cex=0.9,col='darkblue',font=3)
text(x=92, y=240, "TST",cex=0.9,col='darkblue',font=3)
text(x=138, y=240, "HST",cex=0.9,col='darkblue',font=3)
text(x=155, y=240, "FSST",cex=0.9,col='darkblue',font=3)

#12km FSST
h4.12 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_12km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100))

ext4.12 <- cut(h4.12$breaks, c(-Inf,(time_to_strat(1.5,adm_12km,destructive=FALSE)),(time_to_strat(2,adm_12km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h4.12,main= "Last occurrence of taxa at 12 km",                 # histogram with coloring
     sub = "Mass extinction during the Falling Stage system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.12])
  abline(v=(time_to_strat(1.75,adm_12km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=52, "LST",cex=0.9,col='darkblue',font=3)
text(x=13, y=52, "TST",cex=0.9,col='darkblue',font=3)
text(x=22, y=52, "HST",cex=0.9,col='darkblue',font=3)
text(x=50, y=52, "FSST",cex=0.9,col='darkblue',font=3)

#15km FSST
h4.15 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>
  time_to_strat(adm_15km, destructive = FALSE) |>                    # transform into depth domain
  hist(font.axis = 1,                                                # histogram of last occurrences
       font.lab = 3,
       font.sub = 2,
       breaks = seq(from = min_height(adm_15km), to = max_height(adm_15km), length.out = 100))

ext4.15 <- cut(h4.15$breaks, c(-Inf,(time_to_strat(1.5,adm_15km,destructive=FALSE)),(time_to_strat(2,adm_15km,destructive=FALSE)),Inf)) # distinguishing the timing of extinction event
plot(h4.15,main= "Last occurrence of taxa at 15 km",                 # histogram with coloring
     sub = "Mass extinction during the Falling Stage system tract",
     ylab = "Last occurrences per meter",
     xlab = "Stratigraphic height [m]",
     col=c("white","red","white")[ext4.15])
  abline(v=(time_to_strat(1.75,adm_15km,destructive=FALSE)),col="dodgerblue",lwd=3,lty='dashed') # time of peak extinction rate
  abline(v=(time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
  abline(v=(time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
text(x=4, y=54, "LST",cex=0.9,col='darkblue',font=3)
text(x=14, y=54, "TST",cex=0.9,col='darkblue',font=3)
text(x=21.5, y=54, "HST",cex=0.9,col='darkblue',font=3)
text(x=30, y=54, "FSST",cex=0.9,col='darkblue',font=3)
######################################################################

