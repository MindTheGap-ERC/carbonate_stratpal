data_kitten = read.csv("alcap-example_adm_12.csv")
            
#grid: 20km x 50km
#time: 4 Myr
#Sea level from Holland & Patzkowsky 2015

library(admtools)
library(StratPal)
library(ggplot2)
library(dplyr)

#Changes of system tracts
LST_to_TST = 0.5      # Lowstand system tract to Transgressive system tract
TST_to_HST = 1        # Transgressive system tract to Highstand system tract
HST_to_FSST = 1.5     # Highstand system tract to Falling Stage system tract

### Eustatic Sea Level curve
AMP1 = 15                           # Amplitude and period of 3rd order fluctuation
PER1 = 2                            
AMP2 = 2.5                          # Amplitude and period of 4th order fluctuation
PER2 = 0.2
t = seq(1.25, 3.25, by=0.001)       # time frame
t_mod = t[t >= 1.25 & t <= 3.25]    # Subset data
t_mod_shifted = t_mod - min(t_mod)  # shifted from 0 to 2 Myr
sl = (AMP1*sin(2*pi*t/PER1)) + (AMP2*sin(2*pi*t/PER2))    # the curve
plot(t_mod_shifted,sl,type='l',     # Plot curve
     main="Eustatic Sea Level",
     xlab="time [Myr]",
     ylab="[m]")
abline(v=0.5,col="cyan4",lwd=3,lty='dashed')    # separate different system tracts
abline(v=1,col="cyan4",lwd=3,lty='dashed')
abline(v=1.5,col="cyan4",lwd=3,lty='dashed')
text(x=0.2, y=0, "LST",cex=0.9,col='darkblue',font=3)
text(x=0.75, y=12, "TST",cex=0.9,col='darkblue',font=3)
text(x=1.25, y=0, "HST",cex=0.9,col='darkblue',font=3)
text(x=1.75, y=12, "FSST",cex=0.9,col='darkblue',font=3)


# time 4 Myr
t = data_kitten$time..Myr.                              
t_mod = t[t >= 1.25 & t <= 3.25]                          # Subset 1.25 Myr - 3.25 Myr
t_mod_shifted = t_mod - min(t_mod)                        # shift that data to 0 Myr - 2 Myr 

# 20 age depth models (adm). From 1 km to 20 km, every 1 km
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

par(mfrow = c(2, 2))

###  4 km     adm
h4 = data_kitten$adm8..m.      
plot(h4,type='l')                                           # entire runtime (4 Myr)
h4_mod = h4[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h4_mod_shifted = h4_mod - min(h4_mod)                       # shift to 0-2 Myr
adm4_mod = tp_to_adm(t = t_mod, h4_mod)                    
adm_4km = tp_to_adm(t = t_mod_shifted, h4_mod_shifted)     
plot(adm_4km,lwd_acc = 2,lwd_destr = 0)                     # plot 1.25 - 3.25 Myr                              
title(main = "A. ADM 4 km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty="dashed")
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty="dashed")
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty="dashed")
text(x=0.21, y=80, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.75, y=80, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.25, y=40, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.8, y=40, "FSST",cex=0.8,col='darkblue',font=3)


str(adm_4km)
adm_4km_h <- adm_4km$h        # Extract height from adm_4km
adm_4km_t <- adm_4km$t        # and time

df_adm_4km <- data.frame(T_Myr = adm_4km_t,H_4km = adm_4km_h)   # Create a data frame

df_adm_4km <- df_adm_4km %>%                 # Calculate height difference between data points and assign line thickness
  mutate(diff_H_4km = c(NA, diff(H_4km)),
    line_width = ifelse(diff_H_4km == 0, 0.1, 2))

df_segments <- df_adm_4km %>%                # Prepare data for plotting with geom_segment
  mutate(                                    # Each row represents a segment from point i to i+1
    T_Myr_next = lead(T_Myr),
    H_4km_next = lead(H_4km)
  ) %>% filter(!is.na(T_Myr_next))           # Remove the last row (no next point)

# Line plot
ggplot(df_segments) + geom_segment(aes(x=T_Myr, y=H_4km, xend=T_Myr_next,
                                       yend=H_4km_next, size=line_width)) + 
  scale_size_identity() + theme_minimal() + 
  labs(x="Time [Myr]", y="Stratigraphic heigh [m]", 
       title="Age Depth Model at 4 km") + 
  geom_vline(xintercept=LST_to_TST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=TST_to_HST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=HST_to_FSST,color="cyan4",linetype="dashed",size=1) +
  geom_text(aes(x=0.25,y=80,label="LST"),color="darkblue",size=4) + 
  geom_text(aes(x=0.75,y=80,label="TST"),color="darkblue",size=4) +
  geom_text(aes(x=1.25,y=20,label="HST"),color="darkblue",size=4) +
  geom_text(aes(x=1.75,y=20,label="FSST"),color="darkblue",size=4)

#####################################
w4 = water_depth$wd8..m.                                    # Water Depth (wd) at 4 km
plot(w4,type='l',main="complete wd 4km",                    # Entire runtime (4 Myr)
     xlab="time [Myr]",ylab="[m]")
w4_mod = w4[t >= 1.25 & t <= 3.25]                          # modify and shift
w4_mod_shifted = w4_mod - min(w4_mod)
wd4_mod = tp_to_adm(t = t_mod, w4_mod)                    
wd_4km = tp_to_adm(t = t_mod_shifted, w4_mod_shifted) 
plot(wd_4km)                                                # 1.25 - 3.25 Myr
title(main="water depth at 4km",
      ylab="[m]",
      xlab="time [Myr]")
#######################################
###  5 km       adm
h5 = data_kitten$adm10..m. 
plot(h5,type='l')  
h5_mod = h5[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h5_mod_shifted = h5_mod - min(h5_mod)                       # shift to 0-2 Myr
adm5_mod = tp_to_adm(t = t_mod, h5_mod)                    
adm_5km = tp_to_adm(t = t_mod_shifted, h5_mod_shifted)     
plot(adm_5km,lwd_acc = 2,lwd_destr = 0)                     # plot 1.25 - 3.25 Myr                              
title(main = "A. ADM 4 km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)

###  8 km       adm
h8 = data_kitten$adm16..m.
plot(h8,type='l')                                           # Entire runtime (4 Myr)
h8_mod = h8[ t >= 1.25 & t <= 3.25]                         # modify to show 1.25-3.25 Myr
h8_mod_shifted = h8_mod - min(h8_mod)                       # shift to 0-2 Myr
adm8_mod = tp_to_adm(t = t_mod, h8_mod)                     
adm_8km = tp_to_adm(t = t_mod_shifted, h8_mod_shifted)      
plot(adm_8km,lwd_acc = 2,lwd_destr = 0)                     # plot 1.25 - 3.25 Myr
title(main = "B. ADM 8 km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty="dashed")
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty="dashed")
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty="dashed")
text(x=0.22, y=50, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.78, y=50, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.28, y=50, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.8, y=50, "FSST",cex=0.8,col='darkblue',font=3)

###  9 km       adm
h9 = data_kitten$adm18..m.    
plot(h9,type='l')                                           # Entire runtime (4 Myr)
h9_mod = h9[ t >= 1.25 & t <= 3.25]                         # modify to show 1.25-3.25 Myr
h9_mod_shifted = h9_mod - min(h9_mod)                       # shift to 0-2 Myr
adm9_mod = tp_to_adm(t = t_mod, h9_mod)                     
adm_9km = tp_to_adm(t = t_mod_shifted, h9_mod_shifted)      
plot(adm_9km,lwd_acc = 2,lwd_destr = 0)                     # plot 1.25 - 3.25 Myr
title(main = "B. ADM 9 km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty="dashed")
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty="dashed")
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty="dashed")
text(x=0.22, y=50, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.78, y=50, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.28, y=50, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.8, y=50, "FSST",cex=0.8,col='darkblue',font=3)

abline(v=0.52)              # start of slope
abline(v=0.66)              # end
time_to_strat(0.66,adm_9km,destructive=FALSE)-time_to_strat(0.52,adm_9km) # change to depth domain to get slope elevation
time_to_strat(0.52,adm_9km,destructive=FALSE)
adm_9km_h <- adm_9km$h        # Extract height from adm_9km
adm_9km_t <- adm_9km$t        # and time

df_adm_9km <- data.frame(T_Myr = adm_9km_t,H_9km = adm_9km_h)   # Create a data frame

df_adm_9km <- df_adm_9km %>%                 # Calculate height difference between data points and assign line thickness
  mutate(diff_H_9km = c(NA, diff(H_9km)),
         line_width = ifelse(diff_H_9km == 0, 0.1, 2))

df_segments <- df_adm_9km %>%                # Prepare data for plotting with geom_segment
  mutate(                                    # Each row represents a segment from point i to i+1
    T_Myr_next = lead(T_Myr),
    H_9km_next = lead(H_9km)
  ) %>% filter(!is.na(T_Myr_next))           # Remove the last row (no next point)

# Line plot
ggplot(df_segments) + geom_segment(aes(x=T_Myr, y=H_9km, xend=T_Myr_next,
                                       yend=H_9km_next, size=line_width)) + 
  scale_size_identity() + theme_minimal() + 
  labs(x="Time [Myr]", y="Stratigraphic heigh [m]", 
       title="Age Depth Model at 9 km") + 
  geom_vline(xintercept=LST_to_TST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=TST_to_HST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=HST_to_FSST,color="cyan4",linetype="dashed",size=1) +
  geom_text(aes(x=0.25,y=80,label="LST"),color="darkblue",size=4) + 
  geom_text(aes(x=0.75,y=80,label="TST"),color="darkblue",size=4) +
  geom_text(aes(x=1.25,y=20,label="HST"),color="darkblue",size=4) +
  geom_text(aes(x=1.75,y=20,label="FSST"),color="darkblue",size=4)
######################################
w9 = water_depth$wd18..m.                                   # Water depth (wd) at 9 km
plot(w9,type='l',main="complete wd 9km",                    # Entire runtime (4 Myr)
     xlab="time [Myr]",ylab="[m]")
w9_mod = w9[ t >= 1.25 & t <= 3.25]                         # modify and shift
w9_mod_shifted = w9_mod - min(w9_mod)
wd9_mod = tp_to_adm(t = t_mod, w9_mod)                    
wd_9km = tp_to_adm(t = t_mod_shifted, w9_mod_shifted)
plot(wd_9km)                                                # Plot 1.25 - 3.25 Myr
title(main="water depth at 9km",
      ylab="[m]",
      xlab="time [Myr]")
##########################################
### 10km adm
h10 = data_kitten$adm20..m.                           
plot(h10,type='l')                                          # Entire runtime (4 Myr)
h10_mod = h10[ t >= 1.25 & t <= 3.25]
h10_mod_shifted = h10_mod - min(h10_mod)
adm10_mod = tp_to_adm(t = t_mod, h10_mod)                  
adm_10km = tp_to_adm(t = t_mod_shifted, h10_mod_shifted)   
plot(adm_10km,lwd_acc = 2,lwd_destr = 0)                    # 1.25 - 3.25 Myr
title(main = "C. ADM 10 km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')
text(x=0.22, y=80, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.78, y=80, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.29, y=80, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.8, y=80, "FSST",cex=0.8,col='darkblue',font=3)

### 11km adm
h11 = data_kitten$adm22..m.      
plot(h11,type='l')                                          # Entire runtime (4 Myr)
h11_mod = h11[ t >= 1.25 & t <= 3.25]                       # Modify to show 1.25-3.25 Myr
h11_mod_shifted = h11_mod - min(h11_mod)                    # Shift to 0-2 Myr
adm11_mod = tp_to_adm(t = t_mod, h11_mod)                  
adm_11km = tp_to_adm(t = t_mod_shifted, h11_mod_shifted)
plot(adm_11km,lwd_acc = 2,lwd_destr = 0)                    # Plot 1.25 - 3.25
title(main = "C. ADM 11 km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')
text(x=0.22, y=80, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.78, y=80, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.28, y=80, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.9, y=80, "FSST",cex=0.8,col='darkblue',font=3)

abline(v=1.67)            # start of the slope
abline(v=1.78)            # end of the slope
time_to_strat(1.78,adm_11km,destructive=FALSE)-time_to_strat(1.67,adm_11km,destructive=FALSE) # change to depth domain to get slope elevation

adm_11km_h <- adm_11km$h        # Extract height from adm_11km
adm_11km_t <- adm_11km$t        # and time

df_adm_11km <- data.frame(T_Myr = adm_11km_t,H_11km = adm_11km_h)   # Create a data frame

df_adm_11km <- df_adm_11km %>%                 # Calculate height difference between data points and assign line thickness
  mutate(diff_H_11km = c(NA, diff(H_11km)),
         line_width = ifelse(diff_H_11km == 0, 0.1, 2))

df_segments <- df_adm_11km %>%                # Prepare data for plotting with geom_segment
  mutate(                                    # Each row represents a segment from point i to i+1
    T_Myr_next = lead(T_Myr),
    H_11km_next = lead(H_11km)
  ) %>% filter(!is.na(T_Myr_next))           # Remove the last row (no next point)

# Line plot
ggplot(df_segments) + geom_segment(aes(x=T_Myr, y=H_11km, xend=T_Myr_next,
                                       yend=H_11km_next, size=line_width)) + 
  scale_size_identity() + theme_minimal() + 
  labs(x="Time [Myr]", y="Stratigraphic heigh [m]", 
       title="Age Depth Model at 11 km") + 
  geom_vline(xintercept=LST_to_TST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=TST_to_HST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=HST_to_FSST,color="cyan4",linetype="dashed",size=1) +
  geom_text(aes(x=0.25,y=100,label="LST"),color="darkblue",size=4) + 
  geom_text(aes(x=0.75,y=100,label="TST"),color="darkblue",size=4) +
  geom_text(aes(x=1.25,y=100,label="HST"),color="darkblue",size=4) +
  geom_text(aes(x=1.95,y=100,label="FSST"),color="darkblue",size=4)

##################################
w11 = water_depth$wd22..m.                                  # Water depth (wd) at 11 km                                
plot(w11,type='l',main="complete wd 11km",                  # Entire runtime (4 Myr)
     xlab="time [Myr]",ylab="[m]")
w11_mod = w11[ t >= 1.25 & t <= 3.25]                       # Modify and shift
w11_mod_shifted = w11_mod - min(w11_mod)
wd11_mod = tp_to_adm(t = t_mod, w11_mod)                    
wd_11km = tp_to_adm(t = t_mod_shifted, w11_mod_shifted)     
plot(wd_11km)                                               # plot 1.25-3.25 Myr
title(main="water depth at 11km",
      ylab="[m]",
      xlab="time [Myr]")
#####################################
###  12 km    adm
h12 = data_kitten$adm24..m. 
plot(h12,type='l')                                          # Entire runtime (4 Myr)
h12_mod = h12[ t >= 1.25 & t <= 3.25]                       # Modify to show 1.25-3.25 Myr
h12_mod_shifted = h12_mod - min(h12_mod)                    # Sift to 0-2 Myr
adm12_mod = tp_to_adm(t = t_mod, h12_mod)           
adm_12km = tp_to_adm(t = t_mod_shifted, h12_mod_shifted) 
plot(adm_12km,lwd_acc = 2,lwd_destr = 0)                    # 1.25 - 3.25 Myr
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

### 15 km    adm
h15 = data_kitten$adm30..m.   
plot(h15,type='l')                                          # Entire runtime (4 Myr)
h15_mod = h15[ t >= 1.25 & t <= 3.25]                       # Modify to show 1.25-3.25 Myr
h15_mod_shifted = h15_mod - min(h15_mod)                    # Shift to 0-2 Myr
adm15_mod = tp_to_adm(t = t_mod, h15_mod)                  
adm_15km = tp_to_adm(t = t_mod_shifted, h15_mod_shifted)   
plot(adm_15km,lwd_acc = 2,lwd_destr = 0)                    # 1.25 - 3.25 Myr
abline(v=(LST_to_TST),col="cyan4",lwd=3,lty='dashed')
abline(v=(TST_to_HST),col="cyan4",lwd=3,lty='dashed')
abline(v=(HST_to_FSST),col="cyan4",lwd=3,lty='dashed')
title(main = "D. ADM 15 km",
      xlab = "time [Myr]", ylab = "height [m]",
      cex.main = 1.2)
text(x=0.25, y=20, "LST",cex=0.8,col='darkblue',font=3)
text(x=0.75, y=20, "TST",cex=0.8,col='darkblue',font=3)
text(x=1.25, y=10, "HST",cex=0.8,col='darkblue',font=3)
text(x=1.75, y=10, "FSST",cex=0.8,col='darkblue',font=3)

adm_15km_h <- adm_15km$h        # Extract height from adm_15km
adm_15km_t <- adm_15km$t        # and time

df_adm_15km <- data.frame(T_Myr = adm_15km_t,H_15km = adm_15km_h)   # Create a data frame

df_adm_15km <- df_adm_15km %>%                 # Calculate height difference between data points and assign line thickness
  mutate(diff_H_15km = c(NA, diff(H_15km)),
         line_width = ifelse(diff_H_15km == 0, 0.1, 2))

df_segments <- df_adm_15km %>%                # Prepare data for plotting with geom_segment
  mutate(                                     # Each row represents a segment from point i to i+1
    T_Myr_next = lead(T_Myr),
    H_15km_next = lead(H_15km)
  ) %>% filter(!is.na(T_Myr_next))            # Remove the last row (no next point)

# Line plot
ggplot(df_segments) + geom_segment(aes(x=T_Myr, y=H_15km, xend=T_Myr_next,
                                       yend=H_15km_next, size=line_width)) + 
  scale_size_identity() + theme_minimal() + 
  labs(x="Time [Myr]", y="Stratigraphic heigh [m]", 
       title="Age Depth Model at 15 km") + 
  geom_vline(xintercept=LST_to_TST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=TST_to_HST,color="cyan4", linetype="dashed",size=1) + 
  geom_vline(xintercept=HST_to_FSST,color="cyan4",linetype="dashed",size=1) +
  geom_text(aes(x=0.25,y=40,label="LST"),color="darkblue",size=4) + 
  geom_text(aes(x=0.75,y=40,label="TST"),color="darkblue",size=4) +
  geom_text(aes(x=1.25,y=40,label="HST"),color="darkblue",size=4) +
  geom_text(aes(x=1.95,y=40,label="FSST"),color="darkblue",size=4)
###################################
w15 = water_depth$wd30..m.                                    # Water depth (wd) at 15 km                               
plot(w15,type='l',main="complete wd 15km",                    # Entire runtime (4 Myr)
     xlab="time [Myr]")
w15_mod = w15[ t >= 1.25 & t <= 3.25]                         # Modify and shift
w15_mod_shifted = w15_mod - min(w15_mod)
wd15_mod = tp_to_adm(t = t_mod, w15_mod)                    
wd_15km = tp_to_adm(t = t_mod_shifted, w15_mod_shifted)     
plot(wd_15km)                                                 # Plot 1.25 - 3.25
title(main="water depth at 15 km",
      ylab="[m]",
      xlab="time [Myr]")
######################################
###  16 km     adm
h16 = data_kitten$adm32..m.      
plot(h16,type='l')                                           # entire runtime (4 Myr)
h16_mod = h16[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h16_mod_shifted = h16_mod - min(h16_mod)                       # shift to 0-2 Myr
adm16_mod = tp_to_adm(t = t_mod, h16_mod)                    
adm_16km = tp_to_adm(t = t_mod_shifted, h16_mod_shifted)     
plot(adm_16km,lwd_acc = 2,lwd_destr = 0)  

###  18 km     adm
h18 = data_kitten$adm36..m.      
plot(h18,type='l')                                           # entire runtime (4 Myr)
h18_mod = h18[t >= 1.25 & t <= 3.25]                          # modify to show 1.25-3.25 Myr
h18_mod_shifted = h18_mod - min(h18_mod)                       # shift to 0-2 Myr
adm18_mod = tp_to_adm(t = t_mod, h18_mod)                    
adm_18km = tp_to_adm(t = t_mod_shifted, h18_mod_shifted)     
plot(adm_18km,lwd_acc = 2,lwd_destr = 0)  

par(mfrow = c(1, 1))
###############################################################################
install.packages("gridExtra")
library(gridExtra)
### 4 km       constant extinction rate
ext0.4 <- p3(rate = 500, from = min_time(adm_4km), to = max_time(adm_4km)) |>   # constant rate in time domain
  time_to_strat(adm_4km, destructive = FALSE)                                   # transform into depth domain

c_4km <- qplot(ext0.4, geom="histogram",binwidth=1) + coord_flip() + geom_vline(         # plot histogram
  xintercept=time_to_strat(LST_to_TST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "4 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=10,y=100,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=50,y=100,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=85,y=100,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=102,y=100,label="FSST"),color="deepskyblue3",size=4) + 
  theme_minimal()

### 8km        constant extinction rate
ext0.8 <- p3(rate = 500, from = min_time(adm_8km), to = max_time(adm_8km)) |>   # constant rate in time domain
  time_to_strat(adm_9km, destructive = FALSE)                                   # transform into depth domain

qplot(ext0.8, geom="histogram",binwidth=1) + coord_flip() + geom_vline(         # plot histogram
  xintercept=time_to_strat(LST_to_TST,adm_8km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_8km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_8km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "8 km with constant extinction rate",x = "stratigraphic height [m]",y = "Last occurrences") +
  geom_text(aes(x=25,y=80,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=75,y=80,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=120,y=80,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=160,y=80,label="FSST"),color="deepskyblue3",size=4) + 
  theme_minimal()

### 9 km        constant extinction rate
ext0.9 <- p3(rate = 500, from = min_time(adm_9km), to = max_time(adm_9km)) |>   # constant rate in time domain
  time_to_strat(adm_9km, destructive = FALSE)                                   # transform into depth domain

c_9km <- qplot(ext0.9, geom="histogram",binwidth=1) + coord_flip() + geom_vline(         # plot histogram
  xintercept=time_to_strat(LST_to_TST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "9 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=0,y=100,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=70,y=100,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=140,y=100,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=160,y=100,label="FSST"),color="deepskyblue3",size=4) + 
  theme_minimal()

### 10km      constant extinction rate
ext0.10 <- p3(rate = 500, from = min_time(adm_10km), to = max_time(adm_10km)) |>  # constant rate in time domain
  time_to_strat(adm_10km, destructive = FALSE)                                    # transform into depth domain

qplot(ext0.10, geom="histogram",binwidth=1) + coord_flip() + geom_vline(          # plot histogram
  xintercept=time_to_strat(LST_to_TST,adm_10km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_10km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_10km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "10 km with constant extinction rate",x = "stratigraphic height [m]",y = "Last occurrences") +
  geom_text(aes(x=0,y=50,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=50,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=25,y=50,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=100,y=50,label="FSST"),color="deepskyblue3",size=4) + 
  theme_minimal()

### 11 km      constant extinction rate
ext0.11 <- p3(rate = 500, from = min_time(adm_11km), to = max_time(adm_11km)) |>  # constant rate in time domain
  time_to_strat(adm_11km, destructive = FALSE)                                    # transform into depth domain

c_11km <- qplot(ext0.11, geom="histogram",binwidth=1) + coord_flip() + geom_vline(          # plot histogram
  xintercept=time_to_strat(LST_to_TST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "11 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=0,y=52,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=52,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=22,y=52,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=100,y=52,label="FSST"),color="deepskyblue3",size=4) + 
  theme_minimal()

### 12km     constant extinction rate
p3(rate = 500, from = min_time(adm_12km), to = max_time(adm_12km)) |>  # constant rate in time domain
  time_to_strat(adm_12km, destructive = FALSE)                         # transform into depth domain
  #hist(main = "Last occurrence at 12 km with constant extinction rate ",
  #    ylab = "Last occurrence per meter",
  #    xlab = "Stratigraphic height [m]",
  #    breaks = seq(from = min_height(adm_12km), to = max_height(adm_12km), length.out = 100),
  #    col="white")
#abline(v=(time_to_strat(LST_to_TST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(TST_to_HST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#abline(v=(time_to_strat(HST_to_FSST,adm_12km,destructive=FALSE)),col="cyan4",lwd=3,lty='dashed')
#text(x=5, y=30, "LST",cex=0.9,col='darkblue',font=3)
#text(x=15, y=30, "TST",cex=0.9,col='darkblue',font=3)
#text(x=24, y=30, "HST",cex=0.9,col='darkblue',font=3)
#text(x=32, y=30, "FSST",cex=0.9,col='darkblue',font=3)

### 15 km      constant extinction rate
ext0.15 <- p3(rate = 500, from = min_time(adm_15km), to = max_time(adm_15km)) |>  # constant rate in time domain
  time_to_strat(adm_15km, destructive = FALSE)                       # transform into depth domain

qplot(ext0.15, geom="histogram",binwidth=1) + coord_flip() + geom_vline(          # plot histogram
  xintercept=time_to_strat(LST_to_TST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) +
  labs(title = "15 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=4,y=70,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=14,y=70,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=21,y=70,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=32,y=70,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

# 
grid.arrange(c_4km, c_9km, c_11km, c_15km, nrow = 2, ncol = 2)
################################################################################
### Extinction event simulations
# Shape of extinction
LSTsc = c(1,25,1,1,1,1,1,1,1)                                   # LST scenario extinction rate
TSTsc = c(1,1,1,25,1,1,1,1,1)                                   # TST scenario extinction rate
HSTsc = c(1,1,1,1,1,25,1,1,1)                                   # HST scenario extinction rate
FSSTsc = c(1,1,1,1,1,1,1,25,1)                                  # FSST scenario extinction rate
constantrate = c(1,1,1,1,1,1,1,1,1)                             # constant extinction rate
tp = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2)                       # points in time

plot(tp,LSTsc,col="gold1",type='l',lwd=3,                       # plot with all 4 scenarios next to each other
     main="extinction scenarios",
     xlab="time [Myr]",
     ylab="extinction rate [species/m]",
     ylim=c(0,28))
lines(tp, TSTsc, col="gold3",type='l',lwd=3)
lines(tp, HSTsc, col="gold4",type='l',lwd=3)
lines(tp, FSSTsc, col="orange4",type='l',lwd=3)
lines(tp,constantrate,col="black",type='l',lwd=3)
legend(0.2,28,legend=c("LST scen.","TST scen.","HST scen.","FSST scen."),col=c("gold1","gold3","gold4","orange4"),lty=c(1,1,1,1),lwd=3,cex=.85,ncol=4)

### Lowstand Systems Tract extincion - time: 1.25-1.75 Myr (shifted to 0-0.5 Myr)

#4 km LST
ext1.4 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>     # determine position and rate of extinction
  time_to_strat(adm_4km, destructive = FALSE)                           # transform into depth domain

LST_height_4km <- time_to_strat(LST_to_TST,adm_4km,destructive=FALSE)   # maximum height of LST, the end of the extinction event
LST_peak_ext_4km <- time_to_strat((LST_to_TST/2),adm_4km,destructive=FALSE) # height of peak extinction (after 0.25 Myr)
df_ext1.4 <- data.frame(value = ext1.4)                                 # create a data frame

ggplot(df_ext1.4,aes(x=value)) +                                        # plot histogram
  geom_histogram(aes(fill=value <= LST_height_4km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
  xintercept=time_to_strat(LST_to_TST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
    xintercept=time_to_strat(TST_to_HST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
        xintercept=LST_peak_ext_4km,color="red2",linetype="dashed",size=1) +
  labs(title = "LST extinction scenario at 4 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=12,y=145,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=40,y=145,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=80,y=145,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=100,y=145,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#9 km LST
ext1.9 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>   # determine position and rate of extinction
  time_to_strat(adm_9km, destructive = FALSE)                           # transform into depth domain

LST_height_9km <- time_to_strat(LST_to_TST,adm_9km,destructive=FALSE)   #maximum height of LST, the end of the extinction event
LST_peak_ext_9km <- time_to_strat((LST_to_TST/2),adm_9km,destructive=FALSE) #height of peak extinction (after 0.25 Myr)
df_ext1.9 <- data.frame(value = ext1.9)                                 # create a data frame

ggplot(df_ext1.9,aes(x=value)) +                                        # plot histogram
  geom_histogram(aes(fill=value <= LST_height_9km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=LST_peak_ext_9km,color="red2",linetype="dashed",size=1) +
  labs(title = "LST extinction scenario at 9 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=9.7,y=120,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=75,y=120,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=140,y=120,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=155,y=120,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#11km LST
ext1.11 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>    # determine position and rate of extinction
  time_to_strat(adm_11km, destructive = FALSE)                            # transform into depth domain

LST_height_11km <- time_to_strat(LST_to_TST,adm_11km,destructive=FALSE)   #maximum height of LST, the end of the extinction event
LST_peak_ext_11km <- time_to_strat((LST_to_TST/2),adm_11km,destructive=FALSE)  #height of peak extinction (after 0.25 Myr)
df_ext1.11 <- data.frame(value = ext1.11)                                 # create a data frame

ggplot(df_ext1.11,aes(x=value)) +                                         # plot histogram
  geom_histogram(aes(fill=value <= LST_height_11km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=LST_peak_ext_11km,color="red2",linetype="dashed",size=1) +
  labs(title = "LST extinction scenario at 11 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=4,y=130,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=130,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=22,y=130,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=110,y=130,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#15km LST
ext1.15 <- p3_var_rate(x = c(0,0,0.25,0.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>    # determine position and rate of extinction
  time_to_strat(adm_15km, destructive = FALSE)                            # transform into depth domain

LST_height_15km <- time_to_strat(LST_to_TST,adm_15km,destructive=FALSE)   # maximum height of LST, the end of the extinction event
LST_peak_ext_15km <- time_to_strat((LST_to_TST/2),adm_15km,destructive=FALSE)  # height of peak extinction (after 0.25 Myr)
df_ext1.15 <- data.frame(value = ext1.15)                                 # create a data frame

ggplot(df_ext1.15,aes(x=value)) +                                         # plot histogram
  geom_histogram(aes(fill=value <= LST_height_15km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=LST_peak_ext_15km,color="red2",linetype="dashed",size=1) +
  labs(title = "LST extinction scenario at 15 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=2.8,y=150,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=150,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=21,y=150,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=32,y=150,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

### Transgressive Systems Tract extinction - time: 1.75-2.25 Myr (shifted to 0.5-1 Myr)
#4km TST
ext2.4 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>   # determine position and rate of extinction
  time_to_strat(adm_4km, destructive = FALSE)                           # transform into depth domain

TST_height_4km <- time_to_strat(TST_to_HST,adm_4km,destructive=FALSE)   # maximum height of TST, the end of the extinction event
TST_peak_ext_4km <- time_to_strat(((LST_to_TST+TST_to_HST)/2),adm_4km,destructive=FALSE) # height of peak extinction (after 0.75 Myr)
df_ext2.4 <- data.frame(value = ext2.4)                                 # create a data frame

ggplot(df_ext2.4,aes(x=value)) +                                        # plot histogram
  geom_histogram(aes(fill = value >= LST_height_4km & value <= TST_height_4km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=TST_peak_ext_4km,color="red2",linetype="dashed",size=1) +
  labs(title = "TST extinction scenario at 4 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=7,y=50,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=50,y=50,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=85,y=50,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=100,y=50,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#9km TST
ext2.9 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>    # determine position and rate of extinction
  time_to_strat(adm_9km, destructive = FALSE)                           # transform into depth domain

TST_height_9km <- time_to_strat(TST_to_HST,adm_9km,destructive=FALSE)   # maximum height of TST, the end of the extinction event
TST_peak_ext_9km <- time_to_strat(((LST_to_TST+TST_to_HST)/2),adm_9km,destructive=FALSE) # height of peak extinction (after 0.75 Myr)
df_ext2.9 <- data.frame(value = ext2.9)                                 # create a data frame

ggplot(df_ext2.9,aes(x=value)) +                                        # plot histogram
  geom_histogram(aes(fill = value >= LST_height_9km & value <= TST_height_9km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=TST_peak_ext_9km,color="red2",linetype="dashed",size=1) +
  labs(title = "TST extinction scenario at 9 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=6,y=50,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=65,y=50,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=140,y=50,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=155,y=50,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#11km TST
ext2.11 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>    # determine position and rate of extinction
  time_to_strat(adm_11km, destructive = FALSE)                      # transform into depth domain

TST_height_11km <- time_to_strat(TST_to_HST,adm_11km,destructive=FALSE)   # maximum height of TST, the end of the extinction event
TST_peak_ext_11km <- time_to_strat(((LST_to_TST+TST_to_HST)/2),adm_11km,destructive=FALSE) # height of peak extinction (after 0.75 Myr)
df_ext2.11 <- data.frame(value = ext2.11)                                 # create a data frame

ggplot(df_ext2.11,aes(x=value)) +                                        # plot histogram
  geom_histogram(aes(fill = value >= LST_height_11km & value <= TST_height_11km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=TST_peak_ext_11km,color="red2",linetype="dashed",size=1) +
  labs(title = "TST extinction scenario at 11 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=3,y=138,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=138,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=22,y=138,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=100,y=138,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#15km TST
ext2.15 <- p3_var_rate(x = c(0,0.5,0.75,1,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>    # determine position and rate of extinction
  time_to_strat(adm_15km, destructive = FALSE)                       # transform into depth domain

TST_height_15km <- time_to_strat(TST_to_HST,adm_15km,destructive=FALSE)   # maximum height of TST, the end of the extinction event
TST_peak_ext_15km <- time_to_strat(((LST_to_TST+TST_to_HST)/2),adm_15km,destructive=FALSE) # height of peak extinction (after 0.75 Myr)
df_ext2.15 <- data.frame(value = ext2.15)                                 # create a data frame

ggplot(df_ext2.15,aes(x=value)) +                                        # plot histogram
  geom_histogram(aes(fill = value >= LST_height_15km & value <= TST_height_15km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=TST_peak_ext_15km,color="red2",linetype="dashed",size=1) +
  labs(title = "TST extinction scenario at 15 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=5,y=190,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=14,y=190,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=21,y=190,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=30,y=190,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#### Highstand System Tract extinction - time: 2.25-2.75 Myr (shifted to 1-1.5 Myr)
#4km HST
ext3.4 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>   # determine position and rate of extinction
  time_to_strat(adm_4km, destructive = FALSE)                             # transform into depth domain

HST_height_4km <- time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)    # maximum height of HST, the end of the extinction event
HST_peak_ext_4km <- time_to_strat(((TST_to_HST+HST_to_FSST)/2),adm_4km,destructive=FALSE) # height of peak extinction (after 1.25 Myr)
df_ext3.4 <- data.frame(value = ext3.4)                                   # create a data frame

ggplot(df_ext3.4,aes(x=value)) +                                          # plot histogram
  geom_histogram(aes(fill = value >= TST_height_4km & value <= HST_height_4km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=HST_peak_ext_4km,color="red2",linetype="dashed",size=1) +
  labs(title = "HST extinction scenario at 4 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=7,y=150,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=50,y=150,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=80,y=150,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=100,y=150,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#9km HST
ext3.9 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>   # determine position and rate of extinction
  time_to_strat(adm_9km, destructive = FALSE)                             # transform into depth domain
  
HST_height_9km <- time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)    # maximum height of HST, the end of the extinction event
HST_peak_ext_9km <- time_to_strat(((TST_to_HST+HST_to_FSST)/2),adm_9km,destructive=FALSE) # height of peak extinction (after 1.25 Myr)
df_ext3.9 <- data.frame(value = ext3.9)                                   # create a data frame

ggplot(df_ext3.9,aes(x=value)) +                                          # plot histogram
  geom_histogram(aes(fill = value >= TST_height_9km & value <= HST_height_9km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=HST_peak_ext_9km,color="red2",linetype="dashed",size=1) +
  labs(title = "HST extinction scenario at 9 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=6,y=150,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=75,y=150,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=137,y=150,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=155,y=150,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#11km HST
ext3.11 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>    # determine position and rate of extinction
  time_to_strat(adm_11km, destructive = FALSE)                            # transform into depth domain

HST_height_11km <- time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE)    # maximum height of HST, the end of the extinction event
HST_peak_ext_11km <- time_to_strat(((TST_to_HST+HST_to_FSST)/2),adm_11km,destructive=FALSE) # height of peak extinction (after 1.25 Myr)
df_ext3.11 <- data.frame(value = ext3.11)                                   # create a data frame

ggplot(df_ext3.11,aes(x=value)) +                                          # plot histogram
  geom_histogram(aes(fill = value >= TST_height_11km & value <= HST_height_11km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=HST_peak_ext_11km,color="red2",linetype="dashed",size=1) +
  labs(title = "HST extinction scenario at 11 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=4,y=130,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=130,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=22,y=130,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=100,y=130,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#15km HST
ext3.15 <- p3_var_rate(x = c(0,1,1.25,1.5,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>    # determine position and rate of extinction
  time_to_strat(adm_15km, destructive = FALSE)                            # transform into depth domain
  
HST_height_15km <- time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)  # maximum height of HST, the end of the extinction event
HST_peak_ext_15km <- time_to_strat(((TST_to_HST+HST_to_FSST)/2),adm_15km,destructive=FALSE) # height of peak extinction (after 1.25 Myr)
df_ext3.15 <- data.frame(value = ext3.15)                                 # create a data frame

ggplot(df_ext3.15,aes(x=value)) +                                         # plot histogram
  geom_histogram(aes(fill = value >= TST_height_15km & value <= HST_height_15km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=HST_peak_ext_15km,color="red2",linetype="dashed",size=1) +
  labs(title = "HST extinction scenario at 15 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=4,y=120,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=120,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=22,y=120,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=30,y=120,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

### Falling Stage System Tract extinction - time: 2.75-3.25 Myr (shifted to 1.5-2 Myr)
#4km FSST
ext4.4 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>     # determine position and rate of extinction
  time_to_strat(adm_4km, destructive = FALSE)                           # transform into depth domain

HST_height_4km <- time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)  # minimum height of FSST, the onset of the extinction event
FSST_peak_ext_4km <- time_to_strat(((TST_to_HST+TST_to_HST+HST_to_FSST)/2),adm_4km,destructive=FALSE) # height of peak extinction (after 1.75 Myr)
df_ext4.4 <- data.frame(value = ext4.4)                                 # create a data frame

ggplot(df_ext4.4,aes(x=value)) +                                        # plot histogram
  geom_histogram(aes(fill = value >= HST_height_4km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_4km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=FSST_peak_ext_4km,color="red2",linetype="dashed",size=1) +
  labs(title = "FSST extinction scenario at 4 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=7,y=330,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=50,y=330,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=83,y=330,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=105,y=330,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()
time_to_strat(2,adm_4km,destructive=FALSE)-time_to_strat(HST_to_FSST,adm_4km,destructive=FALSE)
time_to_strat(2,adm_15km,destructive=FALSE)-time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)
time_to_strat(1.98,adm_4km,destructive=FALSE)

#9km FSST
ext4.9 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>   # determine position and rate of extinction
  time_to_strat(adm_9km, destructive = FALSE)                           # transform into depth domain

HST_height_9km <- time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE)  # minimum height of FSST, the onset of the extinction event
FSST_peak_ext_9km <- time_to_strat(((TST_to_HST+TST_to_HST+HST_to_FSST)/2),adm_9km,destructive=FALSE) # height of peak extinction (after 1.75 Myr)
df_ext4.9 <- data.frame(value = ext4.9)                                 # create a data frame

ggplot(df_ext4.9,aes(x=value)) +                                        # plot histogram
  geom_histogram(aes(fill = value >= HST_height_9km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_9km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_9km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=FSST_peak_ext_9km,color="red2",linetype="dashed",size=1) +
  labs(title = "FSST extinction scenario at 9 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=5,y=220,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=75,y=220,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=140,y=220,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=158,y=220,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#11km FSST
ext4.11 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>    # determine position and rate of extinction
  time_to_strat(adm_11km, destructive = FALSE)                            # transform into depth domain
  
HST_height_11km <- time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE)  # minimum height of FSST, the onset of the extinction event
FSST_peak_ext_11km <- time_to_strat(((TST_to_HST+TST_to_HST+HST_to_FSST)/2),adm_11km,destructive=FALSE) # height of peak extinction (after 1.75 Myr)
df_ext4.11 <- data.frame(value = ext4.11)                                 # create a data frame

ggplot(df_ext4.11,aes(x=value)) +                                         # plot histogram
  geom_histogram(aes(fill = value >= HST_height_11km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_11km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_11km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=FSST_peak_ext_11km,color="red2",linetype="dashed",size=1) +
  labs(title = "FSST extinction scenario at 11 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=5,y=78,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=78,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=22,y=78,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=90,y=78,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()

#15km FSST
ext4.15 <- p3_var_rate(x = c(0,1.5,1.75,2,2), y = c(1,1,25,1,1), from = 0, to = 2, n = 500, f_max = 50) |>    # determine position and rate of extinction
  time_to_strat(adm_15km, destructive = FALSE)                            # transform into depth domain

HST_height_15km <- time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE)  # minimum height of FSST, the onset of the extinction event
FSST_peak_ext_15km <- time_to_strat(((TST_to_HST+TST_to_HST+HST_to_FSST)/2),adm_15km,destructive=FALSE) # height of peak extinction (after 1.75 Myr)
df_ext4.15 <- data.frame(value = ext4.15)                                 # create a data frame

ggplot(df_ext4.15,aes(x=value)) +                                         # plot histogram
  geom_histogram(aes(fill = value >= HST_height_15km),binwidth=1.5,color="black") + 
  scale_fill_manual(values=c("TRUE"="yellow", "FALSE"="white"), guide="none") + 
  coord_flip() + geom_vline(
    xintercept=time_to_strat(LST_to_TST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
      xintercept=time_to_strat(TST_to_HST,adm_15km,destructive=FALSE),color="cyan4", linetype="dashed",size=1) + geom_vline(
        xintercept=time_to_strat(HST_to_FSST,adm_15km,destructive=FALSE),color="cyan4",linetype="dashed",size=1) + geom_vline(
          xintercept=FSST_peak_ext_15km,color="red2",linetype="dashed",size=1) +
  labs(title = "FSST extinction scenario at 15 km",x = "stratigraphic height [m]",y = "Last occurrences per m") +
  geom_text(aes(x=4,y=150,label="LST"),color="deepskyblue3",size=4) + 
  geom_text(aes(x=13,y=150,label="TST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=22,y=150,label="HST"),color="deepskyblue3",size=4) +
  geom_text(aes(x=28,y=150,label="FSST"),color="deepskyblue3",size=4) +
  theme_minimal()
######################################################################

