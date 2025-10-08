# Import libraries
library(ggplot2)
library(scales)

# Set working directory
setwd("C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/Mutation_rate")

# Define species prefix, either tp or fc - change as needed
PREFIX <- "tp"

# Read in mutation rate table
data <- read.table(paste0("Data/", PREFIX, "_mutRate.txt"), header = FALSE, 
                   col.names = c("Group", "snpRate", "fixRate", "indelRate", "cpRate", "mtRate"))

# Keep sample order as listed in data
ordered_groups <- unique(data$Group)
data$Group <- factor(data$Group, levels = ordered_groups)

# Define sample color palette
if (PREFIX == "tp") {
  colors <- c(
    "T1_9C"  = "#cce5ff",
    "T2_9C"  = "#99ccff",
    "T3_9C"  = "#3399ff",
    "T0_22C" = "#d9f2d9",
    "T1_22C" = "#99e699",
    "T2_22C" = "#33cc33",
    "T1_32C" = "#ffd6d6",
    "T2_32C" = "#ff9999",
    "T3_32C" = "#ff4d4d",
    "T4_32C" = "#cc0000"
  )
} else if (PREFIX == "fc") {
  colors <- c(
    "T1_4C" = "#3399ff",
    "T1_8C" = "#ff4d4d"
  )
} else {
  stop("Unknown PREFIX value. Must be 'tp' or 'fc'.")
}


# Function to plot mutation rate data
plot_box <- function(data, yvar, ylab, colors, filename = NULL, width = 5, height = 5, dpi = 300) {
  p <- ggplot(data, aes(x = Group, y = .data[[yvar]], fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(position = position_jitter(width = 0.2, height=0), size = 2, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scientific_format()) +  # scientific notation
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.1, "cm"),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "",
      y = ylab,
      title = ""
    )
  # Save to file if filename is provided
  if (!is.null(filename)) {
    ggsave(filename, plot = p, width = width, height = height, dpi = dpi, units = "in")
  }
  return(p)
}

# Plot data
plot_box(data, "snpRate", "Novel SNP rate", colors, filename = paste0("Figures/", prefix, "snpRate.jpg"))
plot_box(data, "fixRate", "LOH rate", colors, filename = "Figures/fixRate.jpg")
plot_box(data, "indelRate", "Indel rate", colors, filename = "Figures/indelRate.jpg")
plot_box(data, "cpRate", "cpDNA mutation rate", colors, filename = "Figures/cpRate.jpg") 
plot_box(data, "mtRate", "mtDNA mutation rate", colors, filename = "Figures/mtRate.jpg")
