# Import libraries
library("methylKit")

# Differential methylation analysis
# Very time consuming running the data set en-masse
# so this script allows me to run a group at a time
# or loop through all


# Set working directory
base_dir <- "C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/BS-Seq"
setwd(file.path(base_dir, "bismark_alignments/TP_updated"))
# Set results directory
out_dir <- file.path(base_dir, "differential_methylation/TP_updated")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# =========================
# Function to run differential methylation
# =========================
run_diff_meth <- function(control_files, control_names,
                          treatment_files,  treatment_names,
                          treat_label, assembly="Thaps3", 
                          mincov=5, window_size=50, step_size=50, 
                          diff_threshold=10, qvalue_cutoff=0.01) {
  #message("Running differential methylation for ", treat_label)
  
  # Combine input lists
  file.list <- c(control_files, treatment_files)
  sample.list <- c(control_names, treatment_names)
  treat.groups <- c(rep(0, length(control_files)), rep(1, length(treatment_files)))
  
  # Read and filter data
  myobj <- methRead(file.list,
                    sample.id = sample.list,
                    assembly = assembly,
                    treatment = treat.groups,
                    context = "CpG",
                    mincov = mincov)
  filtered_meth <- filterByCoverage(myobj, lo.count = mincov, hi.count = 100)
  rm(myobj); gc()
  
  # Normalize
  normObj <- normalizeCoverage(filtered_meth, method="median")
  rm(filtered_meth); gc()
  
  # Tile methylation data
  mytiles <- tileMethylCounts(normObj,
                              win.size = window_size,
                              step.size = step_size,
                              cov.bases = 5)
  rm(normObj); gc()
  # Merge
  meth_tiles <- unite(mytiles, destrand = FALSE)
  rm(mytiles); gc()
  
  # Calculate differential methylation
  diff_tiles <- calculateDiffMeth(meth_tiles)
  myDiff.hyper <- getMethylDiff(diff_tiles, difference = diff_threshold, qvalue = qvalue_cutoff, type = "hyper")
  myDiff.hypo  <- getMethylDiff(diff_tiles, difference = diff_threshold, qvalue = qvalue_cutoff, type = "hypo")
  
  # Output results
  write.table(myDiff.hyper, file = file.path(out_dir, paste0(treat_label, "_hyper.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(myDiff.hypo, file = file.path(out_dir, paste0(treat_label, "_hypo.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Finished: ", treat_label)
}


# =========================
# Define sample sets
# =========================
# Define control set
control_files <- list("PRO2022_S1_gDNA_methylKit_format.CpG.txt.gz",
                      "PRO2022_S2_gDNA_methylKit_format.CpG.txt.gz")
control_names <- list("T0_22C_TP1", "T0_22C_TP2")


# Define a named list of treatment sets
treatments <- list(
  "T1_22C" = list(files = list("SAM37194_methylKit_format.CpG.txt.gz", "SAM37195_methylKit_format.CpG.txt.gz",
                               "SAM37196_methylKit_format.CpG.txt.gz", "SAM37197_methylKit_format.CpG.txt.gz",
                               "SAM37198_methylKit_format.CpG.txt.gz"), 
                  names = list("T1_22C_1", "T1_22C_2", "T1_22C_3", "T1_22C_4", "T1_22C_5")),
  "T2_22C" = list(files = list("SAM37209_methylKit_format.CpG.txt.gz", "SAM37210_methylKit_format.CpG.txt.gz", 
                               "SAM37211_methylKit_format.CpG.txt.gz", "SAM37212_methylKit_format.CpG.txt.gz"),
                  names = list("T2_22C_2", "T2_22C_3", "T2_22C_4", "T2_22C_5")),
  "T1_9C"  = list(files = list("SAM37204_methylKit_format.CpG.txt.gz", "SAM37205_methylKit_format.CpG.txt.gz",
                               "SAM37206_methylKit_format.CpG.txt.gz", "SAM37207_methylKit_format.CpG.txt.gz",
                               "SAM37208_methylKit_format.CpG.txt.gz"), 
                  names = list("T1_9C_1", "T1_9C_2", "T1_9C_3", "T1_9C_4", "T1_9C_5")),
  "T2_9C"  = list(files = list("SAM37218_methylKit_format.CpG.txt.gz", "SAM37219_methylKit_format.CpG.txt.gz",
                               "SAM37220_methylKit_format.CpG.txt.gz", "SAM37221_methylKit_format.CpG.txt.gz",
                               "SAM37222_methylKit_format.CpG.txt.gz"), 
                  names = list("T2_9C_1", "T2_9C_2", "T2_9C_3", "T2_9C_4", "T2_9C_5")),
  "T3_9C"  = list(files = list("SAM37223_methylKit_format.CpG.txt.gz", "SAM37224_methylKit_format.CpG.txt.gz",
                               "SAM37225_methylKit_format.CpG.txt.gz", "SAM37226_methylKit_format.CpG.txt.gz",
                               "SAM37227_methylKit_format.CpG.txt.gz"), 
                  names = list("T3_9C_1", "T3_9C_2", "T3_9C_3", "T3_9C_4", "T3_9C_5")),
  "T1_32C" = list(files = list("SAM37199_methylKit_format.CpG.txt.gz", "SAM37200_methylKit_format.CpG.txt.gz",
                               "SAM37201_methylKit_format.CpG.txt.gz", "SAM37202_methylKit_format.CpG.txt.gz",
                               "SAM37203_methylKit_format.CpG.txt.gz"), 
                  names = list("T1_32C_1", "T1_32C_2", "T1_32C_3", "T1_32C_4", "T1_32C_5")),
  "T2_32C" = list(files = list("SAM37213_methylKit_format.CpG.txt.gz", "SAM37214_methylKit_format.CpG.txt.gz",
                               "SAM37215_methylKit_format.CpG.txt.gz", "SAM37216_methylKit_format.CpG.txt.gz",
                               "SAM37217_methylKit_format.CpG.txt.gz"), 
                  names = list("T2_32C_1", "T2_32C_2", "T2_32C_3", "T2_32C_4", "T2_32C_5")),
  "T4_32C" = list(files = list("SAM37229_methylKit_format.CpG.txt.gz", "SAM37230_methylKit_format.CpG.txt.gz",
                               "SAM37232_methylKit_format.CpG.txt.gz"), 
                  names = list("T4_32C_2", "T4_32C_3", "T4_32C_5"))
)


# =========================
# Loop through treatment groups
# =========================
for (label in names(treatments)[1]) {
   treat_set <- treatments[[label]]
   run_diff_meth(control_files, control_names, treat_set$files, treat_set$names, label)
}
