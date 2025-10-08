# Import libraries
library(geosphere)

# Set working directory
setwd("C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/RepairPathways/SOC")

# Read in sample data
samples <- read.table("SOC_coordinates.txt",
                      header=TRUE)

# Generate a vector of all days in the year (non-leap year)
days <- 1:365

# Function to compute median and standard deviation day length for a given latitude
get_daylength_stats <- function(lat) {
  day_lengths <- daylength(lat, days)  # Returns vector of day lengths in hours
  mean_dl <- mean(day_lengths)
  sd_dl <- sd(day_lengths)
  return(c(mean = mean_dl, sd = sd_dl))
}

# Apply function to each location
daylength_stats <- t(apply(samples, 1, function(row) get_daylength_stats(as.numeric(row["Latitude"]))))

# Combine results
result <- cbind(samples, daylength_stats)
#print(result)
                           
# Output results to file
write.table(result, file = "dayLength_by_sample.txt", sep = "\t", quote = FALSE, row.names = FALSE)
