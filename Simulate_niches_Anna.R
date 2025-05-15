library(StratPal) #bio and ecology
library(admtools) #stratigraphy

#task is to create 100 niches described by water depth from CarboKitten model

#Created using chatgpt 

#parameters
set.seed(123)              
n_niches <- 100              # Number of niches
depth_range <- seq(-4.10, 42.44, length.out = 100)  # Water depths from -4.10 to 42.44 meters, min and max across all 4 locations 

# Generate optima more densely toward shallow depth
# Step 1: Create a uniform sequence between 0 and 1
uniform_seq <- seq(0, 1, length.out = n_niches)
# Step 2: Apply a transformation to bias toward 0 (shallow)
# Try squaring or taking sqrt, log, etc.
optima <- -4.10 + (42.44 + 4.10) * (uniform_seq)^2  # Quadratic: more values near minimum depth (for now this means a bias towards -4.10, not sure yet if this should still be around 0)
#other options for bias toward deep: 'sqrt(uniform_seq)' or 'exp(x) / exp(1)'
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
                                     cutoff_val = -4.10)
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

plot(NA, xlim = c(-4.10, 42.44), ylim = c(1, n_niches),
     xlab = "Water depth [m]", ylab = "Niche",
     main = "Depth range of each niche", yaxt = "n")
axis(2, at = 1:n_niches, labels = 1:n_niches, las = 2, cex.axis = 0.5)

for (i in 1:n_niches) {
  if (!is.na(min_depths[i]) && !is.na(max_depths[i])) {
    segments(min_depths[i], i, max_depths[i], i, lwd = 2, col = "blue")
  }
}

#plot(x, numeric_vector(x),
#     type = "l",
#     lwd = 2,
#     xlab = "Water depth [m]",
#     ylab = "Niche optima",
#     main = "Niche optima across different depths")

# Plot to see how optima cluster near shallow depths
#matplot(depth_range, t(niche_matrix[seq(1, 80, by = 8), ]), type = "l",
#        lty = 1, col = rainbow(10), lwd = 2,
#        xlab = "Water Depth (m)", ylab = "Niche Suitability",
#        main = "Niches Clustered at Shallow Depths")
#legend("topright", legend = paste("Niche", seq(1, 80, by = 8)),
#       col = rainbow(10), lty = 1, lwd = 2)

# Plot histogram of niche optima
#hist(optima, breaks = 20, col = "skyblue", border = "white",
#     main = "Distribution of Niche Optima Across Depth",
#     xlab = "Water Depth (m)", ylab = "Number of Niches")