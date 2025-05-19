library(StratPal) #bio and ecology
library(admtools) #stratigraphy

#task is to create 100 niches described by water depth from CarboKitten model

#Created using AI sources and the help of dr. Emilia Jarochowska, Niklas Hohmann, and dr. Charlotte ...

Anna_wd <- read.csv("~/Documents/thesis/example_Anna_wd.csv") #had some problems syncing GitHub and my computer so saved them here locally for now, sorry
Anna_adm <- read.csv("~/Documents/GitHub/stratpal/data/example_Anna_adm.csv")

#as list to try and fix error
Anna_adm_list <- list(
  "t" = Anna_adm$time..Myr.,
  "adm1" = Anna_adm$adm1..m.,
  "adm2" = Anna_adm$adm2..m.,
  "adm3" = Anna_adm$adm3..m.,
  "adm4" = Anna_adm$adm4..m.
)

#parameters
set.seed(124)              
n_niches <- 100              # Number of niches
depth_range <- seq(min(Anna_wd), max(Anna_wd), length.out = 100)  # Water depths from -4.10 to 42.44 meters, min and max across all 4 locations 

# Generate optima more densely toward shallow depth
# Step 1: Create a uniform sequence between 0 and 1
uniform_seq <- seq(0, 1, length.out = n_niches)
# Step 2: Apply a transformation to bias toward 0 (shallow)
# There are optima at wd 0 since those would represent intertidal species, who are quite important when it comes to marine coastal environments (Andrades et al. 2019) and (Bradley et al. 2020)
optima <- max(Anna_wd) * (uniform_seq)^2  # Quadratic: more values near minimum depth (bias towards 0) # -4.10 + (42.44 + 4.10)
#for bias toward shallow: 'log1p(uniform_seq) / log1p(1)'

# Define niche width as a function of depth (e.g., linearly increasing)
min_width <- 2   # Narrowest at shallowest depth
max_width <- 10  # Widest at deepest depth
niche_widths <- min_width + (optima / max(optima)) * (max_width - min_width)

# You can change this to a nonlinear function if you'd like (e.g., log, exp, etc.):
#niche_widths <- min_width + log1p(optima) / log1p(max(optima)) * (max_width - min_width)

# Create list
niches <- list()
for (i in 1:n_niches) {
  niches[[i]] <- StratPal::snd_niche(opt = optima[i],
                                     tol = niche_widths[i],
                                     cutoff_val = 0)
}

niche_val = list(niches)
for (i in 1:100) {
  niche = niches[[i]]
  niche_val[[i]]=niche(depth_range)
} 

niche_mat <- do.call(cbind, niche_val) #converting to matrix (used AI for this line)

#plot niches from list
#used AI to find out how to plot 
matplot(depth_range, niche_mat, type = "l",
        lty = 1, col = rgb(0,0,1,0.2),
        xlab = "Water depth [m]",
        ylab = "Niche value",
        main = "Niche optima across different depths")

#(using AI) plotting them as depth ranges to show more clearly the range and abundance differences across the depths
#note: the ranges continue after 80 for the highest numbered niches, because their optima fall close to or around 80, which means although it does not show they still have the broadest suitable depth range 
threshold <- 0.005  # Set your threshold; chose 0.005 since that means there is only a 0.5% likelihood of finding the organism at that depth

n_niches <- length(niches)
min_depths <- numeric(n_niches)
max_depths <- numeric(n_niches)

for (i in 1:n_niches) {
  values <- niches[[i]](depth_range)
  # Use a threshold to define "present" (e.g., >0 or >0.01)
  present <- which(values > threshold)
  if (length(present) > 0) {
    min_depths[i] <- depth_range[min(present)]
    max_depths[i] <- depth_range[max(present)]
  } else {
    min_depths[i] <- NA
    max_depths[i] <- NA
  }
}

plot(NA, xlim = c(min(Anna_wd), max(Anna_wd)), ylim = c(1, n_niches),
     xlab = "Water depth [m]", ylab = "Niche",
     main = "Depth range of each niche", yaxt = "n")
axis(2, at = 1:n_niches, labels = 1:n_niches, las = 2, cex.axis = 0.5)

for (i in 1:n_niches) {
  if (!is.na(min_depths[i]) && !is.na(max_depths[i])) {
    segments(min_depths[i], i, max_depths[i], i, lwd = 2, col = "blue")
  }
}

# Plot histogram of niche optima
hist(optima, breaks = 20, col = "blue", border = "white",
     main = "Distribution of Niche Optima Across Depth",
     xlab = "Water Depth (m)", ylab = "Number of Niches")

#adapted from:
vignette("event_data")

#change numbers to 1, 2, 3, or 4 depending on which location
#Niche modeling
#how gradient changes with time (in this case water depth)
t = Anna_wd$time..Myr.           # time steps of the model
wd = Anna_wd$wd1..m.   # water depth 2 km offshore at model time steps
gc = approxfun(t, wd)         # define function that defines how the gradient changes with time (gc = *G*radient *C*hange)

plot(t, gc(t), 
     type = "l", 
     xlab = "Time [Myr]",
     ylab = "Water depth [m]",
     main = "Water depth 1.5 km offshore")

#niche model in time dimension

set.seed(124)

niches_applied_t <- list()
for (i in 1:100) {
  niches_applied_t[[i]] <- p3(rate = 300, from = min(t), to = max(t)) |> 
    apply_niche(niche_def = niches[[i]], gc = gc)
}

#histograms
hist_list_t <- lapply(niches_applied_t, function(x) {
  x_num <- as.numeric(x)
  x_num <- x_num[!is.na(x_num)]
  if (length(x_num) > 1 && length(unique(x_num)) > 1) {
    hist(x_num, plot = FALSE)
  } else {
    NULL
  }
})

# Remove NULLs
hist_list_t <- Filter(Negate(is.null), hist_list_t)

# Combine into a data frame
hist_data_t <- do.call(rbind, lapply(seq_along(hist_list_t), function(i) {
  h <- hist_list_t[[i]]
  data.frame(
    mid = h$mids,
    count = h$counts,
    niche = as.factor(i)
  )
}))

library(ggplot2)
ggplot(hist_data_t, aes(x = mid, y = count, fill = niche)) +
  geom_density(aes(y = after_stat(count)), alpha = 0.15) +
  labs (x = "Time [Myr]", y = "#Fossils")

#niche model in terms of stratigraphy 

niches_applied_strat <- list()
for (i in 1:3) {
  niches_applied_strat[[i]] <- p3(rate = 300, from = min(t), to = max(t)) |> 
    apply_niche(niche_def = niches[[i]], gc = gc) |>                    
    time_to_strat(Anna_adm_list$adm1, t, destructive = TRUE, out_dom_val_h = "default") |>                     # transform into strat. domain, destroy fossils that coincide with hiatuses 
    hist(xlab = "Stratigraphic height [m]",                       
         main = "Fossil abundance 1.5 km from shore",
         ylab = "# Fossils",
         breaks = seq(from = min(Anna_adm), to = max(Anna_adm), length.out = 60))
}

#water depth as presented by the stratigraphy  
list("t" = t, "y" = wd) |>    # create list with time - water depth information
  time_to_strat(Anna_adm_list$adm1) |>   # transform into the strat. domain
  plot(orientation = "lr",    # plot water depth information in the stratigraphic domain
       type = "l",
       xlab = "Stratigraphic position [m]",
       ylab = "Water depth [m]",
       main = "Water depth in section 1.5 km from shore")

