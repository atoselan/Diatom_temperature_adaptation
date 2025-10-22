# Import libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)

# Set working directory
setwd("C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/Coverage/CNV_geneCentric/Per_sample_cnv")

# Min copy number change cutoff
min_cnv <- 0.25

# Read in cnv data
myDat <- read.table("FC_patched_over100k_contigs_normCov_perGene.txt",
                    sep="\t",
                    header=TRUE)
# Initialise column for relative position
myDat$RelPos <- 0

# Define subset of scaffolds to keep
n75_scaffolds <- c("scaffold_1a", "scaffold_1b", "scaffold_2", "scaffold_3", "scaffold_4",
                   "scaffold_5", "scaffold_6", "scaffold_7", "scaffold_8", "scaffold_9",
                   "scaffold_10", "scaffold_11", "scaffold_12", "scaffold_13", "scaffold_14",
                   "scaffold_15", "scaffold_16", "scaffold_17", "scaffold_18", "scaffold_19",
                   "scaffold_20", "scaffold_21", "scaffold_22", "scaffold_23", "scaffold_24",
                   "scaffold_25")

# Update scaffold_1. Split into scaffold_1a and scaffold_1b
# around the position of the predicted telomere motif

# Function to update scaffold names and positions
update_scaffold_data <- function(df) {
  df %>%
    mutate(
      scaffold = case_when(
        scaffold == "scaffold_1" & start < split_position ~ "scaffold_1a",
        scaffold == "scaffold_1" & start >= split_position ~ "scaffold_1b",
        TRUE ~ scaffold
      ),
      start = ifelse(scaffold == "scaffold_1b", start - (split_position - 1), start),
      end = ifelse(scaffold == "scaffold_1b", end - (split_position - 1), end)
    )
}
# Define split position
split_position <- 7699806  

# Apply the function to all datasets
myDat_updated <- update_scaffold_data(myDat)

# Find the chromosome length for each scaffold
# Get the last row for each chromosome
scaffold_ends <- myDat_updated %>%
  group_by(scaffold) %>%
  summarise(max_end = max(end))
scaffold_ends <- as.data.frame(scaffold_ends)

# Add relative position to data
for (row in 1:nrow(myDat_updated))
{
  c <- myDat_updated$scaffold[row]
  c_length <- scaffold_ends[scaffold_ends$scaffold == c, 2]
  start <- myDat_updated$start[row]
  stop  <- myDat_updated$end[row]
  # Calculate relative position as middle of gene
  myDat_updated$RelPos[row] <- (((start + stop) /2) / c_length) * 100
}

# Add relative position percentile group
myDat_updated$RelGroup <- floor(myDat_updated$RelPos / 10) + 1

# Filter data for just N75 scaffolds
myDat_N75 <- myDat_updated[myDat_updated$scaffold %in% n75_scaffolds, ]



###
### cnv scatterplot
###

# Subset data by temperature and timepoint
ctrl <- subset(myDat_N75, sample=="FC_T0_4C_S1" | sample=="FC_T0_4C_S2" | sample=="FC_T0_4C_S3")
low_tmp <- subset(myDat_N75, sample=="FC_T66_4C_R1" | sample=="FC_T66_4C_R2" | sample=="FC_T66_4C_R3")
high_tmp <- subset(myDat_N75, sample=="FC_T47_8C_R1" | sample=="FC_T47_8C_R2" | sample=="FC_T47_8C_R3")

# Calculate median copy number
ctrl_median <- ctrl %>%
  group_by(scaffold, start, end, RelPos) %>%
  summarise(median_norm_cov = median(norm_cov, na.rm = TRUE), .groups = "drop")
low_median <- low_tmp %>%
  group_by(scaffold, start, end, RelPos) %>%
  summarise(median_norm_cov = median(norm_cov, na.rm = TRUE), .groups = "drop")
high_median <- high_tmp %>%
  group_by(scaffold, start, end, RelPos) %>%
  summarise(median_norm_cov = median(norm_cov, na.rm = TRUE), .groups = "drop")

# Add control median to each subset
low_median$Ctrl_Median <- ctrl_median$median_norm_cov
high_median$Ctrl_Median <- ctrl_median$median_norm_cov
# Create T1 high temp vs T1 low temp dataset
high_vs_low_median <- high_median
high_vs_low_median$Ctrl_Median <- low_median$median_norm_cov

# Change column name to match plot function
colnames(low_median)[colnames(low_median) == "median_norm_cov"] <- "norm_cov"
colnames(high_median)[colnames(high_median) == "median_norm_cov"] <- "norm_cov"
colnames(high_vs_low_median)[colnames(high_vs_low_median) == "median_norm_cov"] <- "norm_cov"

# Calculate difference
low_median$Diff  <- low_median$norm_cov - low_median$Ctrl_Median
high_median$Diff  <- high_median$norm_cov - high_median$Ctrl_Median
high_vs_low_median$Diff <- high_vs_low_median$norm_cov - high_vs_low_median$Ctrl_Median


# Function to scatterplot copy number compared to T0
plotCnv <- function(dataSet, sampleName, myCol, ylabel="", tempLabel=""){
  
  plot(dataSet$Ctrl_Median~dataSet$norm_cov,
       xlim=c(0, 6.5),
       ylim=c(0, 6.5),
       axes=FALSE,
       ann=FALSE,
       pch=16,
       col=alpha(myCol, 0.4),
       cex=1.25, cex.lab=1.5)
  axis(1, at = seq(0, 6.5, by = 0.5), labels = seq(0, 6.5, by = 0.5), cex.axis = 1.5)
  axis(2, at = seq(0, 6.5, by = 0.5), labels = seq(0, 6.5, by = 0.5), cex.axis = 1.5)
  # Add x-axis label
  mtext(sampleName, side=1, line=3, cex=1)
  # Add y-axis label
  mtext(ylabel, side=2, line=3, cex=1)
  # Add temperature label
  mtext(tempLabel, side=4, line=1, cex=1)
  
  abline(a = 0, b = 1,
         lty="dashed",
         lwd=1,
         col="black")
  
  # Plot CNV loss - TN - T0 <= -0.25
  points(subset(dataSet, Diff <= -min_cnv)$Ctrl_Median~subset(dataSet, Diff <= (min_cnv*-1))$norm_cov,
         col=alpha("#fdae61", 0.5), pch=16, cex=1.25)
  points(subset(dataSet, Diff <= -min_cnv)$Ctrl_Median~subset(dataSet, Diff <= (min_cnv*-1))$norm_cov,
         col=alpha("#d7191c", 0.5), pch=21, cex=1.25)
  
  # Plot CNV gain - TN - T0 >= 0.25
  points(subset(dataSet, Diff >= min_cnv)$Ctrl_Median~subset(dataSet, Diff >= min_cnv)$norm_cov,
         col=alpha("#abd9e9", 0.5), pch=16, cex=1.25)
  points(subset(dataSet, Diff >= min_cnv)$Ctrl_Median~subset(dataSet, Diff >= min_cnv)$norm_cov,
         col=alpha("#2c7bb6", 0.5), pch=21, cex=1.25)
}

# Print plot to file
pdf("C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/Coverage/Summary_figures/FC_cnv_geneCentric_plot_cnv025.pdf",
    width=15,height=6)
par(mfrow=c(1,3))
plotCnv(low_median[low_median$Ctrl_Median>0,], "T1_4C", "#cccccc", ylabel="T0_4C")
plotCnv(high_median[high_median$Ctrl_Median>0,], "T1_8C", "#cccccc", ylabel="T0_4C")
plotCnv(high_vs_low_median[high_vs_low_median$Ctrl_Median>0,], "T1_8C", "#cccccc", ylabel="T1_4C")
dev.off()



###
### Bubble plot of cnv losses and gains
###

# Summarize gains and losses per scaffold
bubble_data1 <- low_median %>%
  mutate(Type = case_when(
    Diff >= 0.25 ~ "Gain",
    Diff <= -0.25 ~ "Loss",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Type)) %>%
  dplyr::count(scaffold, Type)  # Count number of gains/losses per scaffold
# Ensure scaffold order remains the same as the n75 list
bubble_data1 <- bubble_data1 %>%
  mutate(scaffold = factor(scaffold, levels = n75_scaffolds))

bubble_data2 <- high_median %>%
  mutate(Type = case_when(
    Diff >= 0.25 ~ "Gain",
    Diff <= -0.25 ~ "Loss",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Type)) %>%
  dplyr::count(scaffold, Type)  # Count number of gains/losses per scaffold
bubble_data2 <- bubble_data2 %>%
  mutate(scaffold = factor(scaffold, levels = n75_scaffolds))


# Subset into losses or gains - just change the type filter
gain_data1 <- bubble_data1 %>% filter(Type == "Gain")
gain_data2 <- bubble_data2 %>% filter(Type == "Gain")

# Add a new column to identify which dataset each row came from
gain_data1 <- gain_data1 %>% mutate(Source = "Low")
gain_data2 <- gain_data2 %>% mutate(Source = "High")

# Combine them
combined_gain_data <- bind_rows(gain_data1, gain_data2)

# Abbreviate scaffold names
# Remove prefix and convert to factor in the desired order
combined_gain_data <- combined_gain_data %>%
  mutate(
    scaffold = str_remove(as.character(scaffold), "^scaffold_"),
    scaffold = factor(scaffold, levels = str_remove(n75_scaffolds, "^scaffold_"))
  )

# Define colors and sample labels
colors <- c(rep(alpha("#abd9e9",0.5),1)) # light blue for gains
#colors <- c(rep(alpha("#fdae61",0.5),1)) # orange colour for losses
sampleLabels <- c("T1_4C", "T1_8C")

# Create plot
p <- ggplot(combined_gain_data, aes(x=Source, y=scaffold)) + 
  geom_point(aes(size=n, fill=Type), alpha=0.75, shape=21) +
  scale_size_continuous(limits = c(0, 750), range = c(0,15), breaks = c(1,50,100,250,750)) +
  geom_vline(xintercept=1.5, size=1, color="black", linetype="dotted") +
  #geom_vline(xintercept=2, size=1, color="black", linetype="dotted") +
  labs(x="", y="Scaffold", size="Frequency", fill="") +
  scale_x_discrete(labels=sampleLabels, position="top") +
  theme(legend.key=element_blank(),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", vjust = 0.5, hjust = 0.5), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") + 
  scale_fill_manual(values = colors , guide = FALSE)

# Save plot to pdf
ggsave(
  "C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/Coverage/CNV_geneCentric/FC_highVlow_perChrom_gains_025.png",
  width = 4,
  height = 8,
  dpi = 300,       # high quality
  units = "in"     # width and height in inches
)


###
### Plot ratio of gains to losses over relative position
###

plot_cnv_ratio <- function(data, min_cnv = 0.25, 
                           filename = NULL,
                           width = 6, height =5, res = 300,
                           color = "black",
                           main_title = "CNV ratio plot") {
  
  # If a filename is given, open a PNG device
  if (!is.null(filename)) {
    png(filename, width = width, height = height, units = "in", res = res)
    on.exit(dev.off())  # Ensure device closes when function exits
  }
  
  par(mfrow = c(1, 1))  # Single plot per call
  
  # Define bins (0–100%)
  bins <- seq(0, 100, by = 10)
  bin_midpoints <- bins[-length(bins)] + diff(bins) / 2
  
  # Extract positions for gains and losses
  gains <- data[data$Diff >= min_cnv, ]$RelPos
  losses <- data[data$Diff <= -min_cnv, ]$RelPos
  
  # Bin data
  gain_bins <- cut(gains, breaks = bins, include.lowest = TRUE, right = FALSE)
  loss_bins <- cut(losses, breaks = bins, include.lowest = TRUE, right = FALSE)
  
  # Count per bin
  gain_counts <- as.numeric(table(gain_bins))
  loss_counts <- as.numeric(table(loss_bins))
  
  # Calculate log2 ratio (handle divide-by-zero with NA)
  ratio <- log2(ifelse(loss_counts == 0, NA, gain_counts / loss_counts))
  
  # Build dataframe for reference (optional to return)
  ratio_df <- data.frame(
    bin = levels(gain_bins),
    midpoint = bin_midpoints,
    gains = gain_counts,
    losses = loss_counts,
    ratio = ratio
  )
  
  # Plot ratio vs. relative position
  plot(ratio_df$ratio, type = "b",
       main = main_title,
       cex.main = 1.5,
       ylab = "log2(gains/losses)",
       xlab = "Relative scaffold\nposition",
       cex.lab = 1.2,
       ylim = c(-2, 2),
       xaxt = "n",
       yaxt = "n",
       col = color,
       pch = 16,
       cex = 2.5)
  # Add axes and reference line
  axis(side = 1, at = c(1, 5.5, 10), 
       labels = c(0, 0.5, 1), 
       cex.axis = 1.2)
  axis(side = 2, at = c(-2, -1, 0, 1, 2), 
       labels = c("-2.0", "-1.0", "0.0", "1.0", "2.0"),
       cex.axis = 1.2)
  abline(h = 0, lty = 2)
  
  # Return ratio dataframe for further use if needed
  return(ratio_df)
}

# Plot low temperature T1 vs T0
result <- plot_cnv_ratio(
  data = low_median,
  min_cnv = 0.2,
  color = "#1f78b4",
  main_title = "T1_4C vs T0_4C"
)
# Plot high temperature T1 vs T0
result <- plot_cnv_ratio(
  data = high_median,
  min_cnv = 0.2,
  color = "#e31a1c",
  main_title = "T1_8C vs T0_4C"
)
# Plot high vs low T1s
result <- plot_cnv_ratio(
  data = high_vs_low_median,
  min_cnv = 0.2,
  color = "#ff7f00",
  main_title = "T1_8C vs T1_4C"
)
