# Import libraries
library(pheatmap)
library(mgcv)

# Set data directory
hr_og_dir <- "C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/BRCA2/Environmental_BRCA2/HDR_orthogroups"

# Read in data
myDat <- read.table(paste(hr_og_dir, "/", "HDR_OG_updated_normalised_counts.txt", sep=""),
                    header=TRUE,
                    sep="\t")
# Filter out longitudinal transect samples
myDat <- subset(myDat, latitude<=78.62)
# Remove ortholog group part from column names from the 4th column onward
colnames(myDat)[10:ncol(myDat)] <- sub(".*\\.", "", colnames(myDat)[10:ncol(myDat)])
# Set sample IDs as rownames
rownames(myDat) <- myDat[[1]]

# Create colour palette for main heatmap
my_colors <- colorRampPalette(c("white", "#000000"))(10)

# Set up annotation data frame
annotation_row <- data.frame(
  Latitude = myDat$latitude,
  DayLength_SD = myDat$sd_dayLength
)
# Ensure annotation row names match data row names
rownames(annotation_row) <- rownames(myDat)


# Bin latitude: from -80 to 80 in 20° intervals
annotation_row$Latitude <- cut(
  myDat$latitude,
  breaks = seq(-80, 80, by = 20),
  include.lowest = TRUE,
  right = FALSE,  # means intervals are like [-80, -60)
  labels = c("60°N to 80°N", "40°N to 60°N", "20°N to 40°N", "0° to 20°N",
             "0° to 20°S", "20°S to 40°S","40°S to 60°S","60°S to 80°S")
)
# Bin day length standard deviation from 0 to 12 hours by 2 hours
annotation_row$DayLength_SD <- cut(
  myDat$sd_dayLength,
  breaks = seq(0, 12, by = 2),
  include.lowest = TRUE,
  right = FALSE,
  labels = paste0(seq(0, 10, by = 2), "–", seq(2, 12, by = 2))
)

# Assign colors to annotation variables
lat_colors <- setNames(
  colorRampPalette(c("#313695", "#74add1", "#abd9e9",
                     "white", "white",
                     "#abd9e9", "#74add1", "#313695"))(length(levels(annotation_row$Latitude))),
  levels(annotation_row$Latitude))
day_colors <- setNames(
  colorRampPalette(c("white", "#fee08b", "#fdae61", "#f46d43", "#d53e4f"))(length(levels(annotation_row$DayLength_SD))),
  levels(annotation_row$DayLength_SD)
)
# Combine
annotation_colors <- list(
  Latitude = lat_colors,
  DayLength_SD = day_colors
)
  
# Plot heatmap
png("C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/BRCA2/Environmental_BRCA2/HDR_orthogroups/HR_SOC_heatmap_2annot.png", 
    width = 2000, height = 3000, res = 300)
  # Display heatmap
  pheatmap(
    log(myDat[-c(1:11)]+1),
    #breaks = my_breaks,
    color = my_colors,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = TRUE,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    fontsize_row = 4,
    fontsize_col = 10,
    border_color = "black",
    main = "HDR"
  )
dev.off()

  
### GAM analysis vs Latitude

row_means <- rowMeans(log(myDat[-c(1:11)]+1))
latitudes <- myDat$latitude
  
# Fit GAM model
gam_fit <- gam(row_means ~ s(latitudes))
# Get summary
gam_sum <- summary(gam_fit)

# Extract useful parameters
N <- gam_sum$n             # number of observations
adj_r2 <- round(gam_sum$r.sq, 2)   # adjusted R²
edf <- round(gam_sum$s.table[1, "edf"], 3)  # EDF for the smooth term
p_val <- gam_sum$s.table[1, "p-value"]      # p-value for the smooth term
p_text <- ifelse(p_val < 2e-16, "<2e-16", formatC(p_val, format = "e", digits = 2))
  
# Predict over a sequence of latitude values for smooth curve
lat_seq <- seq(min(latitudes), max(latitudes), length.out = 200)
pred <- predict(gam_fit, newdata = data.frame(latitudes = lat_seq), se.fit = TRUE)

# Plot the data
png("C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/BRCA2/Environmental_BRCA2/HDR_orthogroups/HDR_SOC_latitude_GAM.png", 
    width = 1500, height = 1500, res = 300)
  plot(latitudes, row_means, 
       pch = 16,
       cex = 2,
       cex.lab = 1.2,
       cex.axis = 1.2,
       col = "#f46d43",
       xlab = "Latitude", ylab = "Mean log(TPM)",
       main = "Mean HDR vs Latitude")
  points(latitudes, row_means, 
         pch = 21, 
         cex=2, 
         col = "black")
  # Create legend text
  legend_text <- c(
    "GAM (Gaussian)",
    paste0("N = ", N),
    paste0("Adj. R² = ", adj_r2),
    paste0("P = ", p_text),
    paste0("EDF = ", edf)
  )
  # Add text box
  legend("top", legend = legend_text, bty = "n", cex = 1)
  # Add the GAM fit line
  lines(lat_seq, pred$fit, col = "blue", lwd = 2)
  # Add confidence interval
  lines(lat_seq, pred$fit + 2 * pred$se.fit, col = "blue", lty = 2)
  lines(lat_seq, pred$fit - 2 * pred$se.fit, col = "blue", lty = 2)
dev.off()
