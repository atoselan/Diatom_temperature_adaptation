###
### Normalise matches to orthogroups from SOC metatranscriptomes to TPM
###


# Define data locations
hmm_folder <- "C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/BRCA2/Environmental_BRCA2/HDR_orthogroups/HDR_SOC_hmm_results"
out_folder <- "C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/BRCA2/Environmental_BRCA2/HDR_orthogroups/HDR_SOC_hmm_normalised"
map_folder <- "C://Users/andre/OneDrive - University of East Anglia/Sea_of_change_data/metat/self_map"
bitscore_file <- "C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/BRCA2/Environmental_BRCA2/HDR_orthogroups/bitscore_thresholds.csv"

# Load table of orthogroup-specific bitscore thresholds
# Should have columns: OG, bitscore
bitscore_table <- read.table(bitscore_file, header=TRUE, stringsAsFactors = FALSE)

# Set fairly weak evalue cutoff, just to remove weak hits
# Main filtering will be done by bitscore
eval_cutoff <- 1e-5

# Get list of all orthogroups (i.e., subdirectories)
OG_dirs <- list.dirs(hmm_folder, recursive = FALSE, full.names = FALSE)


# File and sample prefixes
sample_prefix <- c("ARK-5M-1", "ARK-5M-2", "ARK-5M-3", "ARK-7M-1", "ARK-7M-2",
                   "ARK-7M-3", "ANT-15_ARK-20-1", "ANT-15_ARK-20-2", "ANT-1", "ANT-10",
                   "ANT-11", "ANT-12", "ANT-13", "ANT-14", "ANT-16", "ANT-17",
                   "ANT-18", "ANT-19", "ANT-2", "ANT-3", "ANT-4", "ANT-5",
                   "ANT-8", "ARK-118M", "ARK-122M", "ARK-130M", "ARK-132M",
                   "ARK-132S", "ARK-135M", "ARK-149", "ARK-15M", "ARK-1M",
                   "ARK-1S", "ARK-20M", "ARK-27M", "ARK-2M", "ARK-2S",
                   "ARK-37M", "ARK-3M", "ARK-3S", "ARK-51M", "ARK-55S",
                   "ARK-61M", "ARK-63M", "ARK-87M", "ARK-8M", "ANT-15_ARK-20-2",
                   "PS103-1.R1", "PS103-1.R2", "PS103-1.R3", "PS103-11",
                   "PS103-14", "PS103-17", "PS103-2.R1", "PS103-2.R2",
                   "PS103-2.R3", "PS103-23", "PS103-29", "PS103-3.R1",
                   "PS103-3.R2", "PS103-3.R3", "PS103-34", "PS103-4.R1",
                   "PS103-4.R2", "PS103-4.R3", "PS103-5.R1", "PS103-5.R2",
                   "PS103-5.R3", "PS103-59", "PS103-6", "PS103-7", "S0_C1_B9",
                   "S1_C1_B19", "S11_C1_B13", "S13_C1_B13", "S15_C1_B20",
                   "S17_C1_B21", "S18_C1_B14", "S2_C1_B8", "S2_C18_B9",
                   "S20_C1_B20", "S21_C1_B12", "S23_C1_B24", "S25_C1_B22",
                   "S27_C1_B23", "S29_C1_B24", "S3_C1_B8", "S30_C1_B23",
                   "S32_C1_B24", "S5_C1_B9", "S5_C27_B21", "S7_C1_B21",
                   "S9_C1_B7")
hmm_prefix <- c("164863", "164862", "164799", "164797", "164798",
                "164791", "164794", "164796", "164845", "164819",
                "164848", "164826", "164800", "164827", "164828",
                "164801", "164829", "164802", "164847", "164846",
                "164821", "164817", "164850", "164861", "164841",
                "164807", "164842", "164843", "164820", "164844",
                "164814", "164831", "164832", "164811", "164836",
                "164833", "164834", "164815","164806", "164805",
                "164837", "164838", "164816", "164839", "164840",
                "164835", "164795", "203299", "203307", "203306",
                "200440", "200443", "200441", "203305", "203302",
                "203303", "200438", "200442", "203301", "202327",
                "203178", "200439", "203236", "203298", "202539",
                "203083", "203235", "203300", "203288", "203304",
                "203234", "164851", "164852", "164808", "164810",
                "164849", "164830", "164813", "164853", "164854",
                "164855", "164856", "164809", "164857", "164860",
                "164858", "164824", "164822", "164823", "164812",
                "164825", "164803", "164804")
selfmap_prefix <- c("205070", "204834", "204857", "204847", "204848",
                    "204744", "205093", "205073", "204742", "204861",
                    "204719", "204844", "204843", "205090", "205084",
                    "204856", "205091", "205071", "204723", "204730",
                    "204838", "205076", "204743", "205086", "205078",
                    "204862", "204715", "204741", "205074", "204729",
                    "204873", "204869", "204731", "204840", "204747",
                    "204718", "204725", "204854", "204863", "204850",
                    "204837", "204870", "204839", "204877", "205080",
                    "204713", "205094", "203775", "203857", "203863",
                    "202402", "202397", "202422", "203867", "203778",
                    "203777", "202419", "202389", "203864", "202523",
                    "203858", "202394", "210577", "203866", "203612",
                    "203611", "203862", "203860", "203861", "203774",
                    "203859", "204702", "204845", "204867", "204714",
                    "204749", "204880", "204864", "204878", "204835",
                    "204849", "204842", "205088", "204846", "204879",
                    "204853", "204855", "204859", "204874", "205077",
                    "204866", "205083", "204881")
# Combine sample names and file prefixes into data frame
details_tab <- data.frame(sample_prefix,
                          hmm_prefix,
                          selfmap_prefix,
                          stringsAsFactors = FALSE)


# Loop through ortholog groups and normalise
for (OG in OG_dirs) {
  
  # Get bitscore threshold for this OG
  bitscore_thresh <- bitscore_table$Bitscore_10pct[bitscore_table$Orthogroup == OG]
  if (length(bitscore_thresh) == 0) {
    message("### No bitscore threshold found for ", OG, ", skipping.")
    next
  }
  
  # Initialise list to store per-sample TPM sums (default to 0)
  result_list <- setNames(rep(0, nrow(details_tab)), details_tab$sample_prefix)
  
  # Loop through samples
  for (i in 1:nrow(details_tab)) {
    # Get sample info
    sample_name <- details_tab[i, 'sample_prefix']
    hmm_file <- file.path(hmm_folder, OG, paste0(details_tab[i, 'hmm_prefix'], "_", OG, ".domtblout.gz"))
    map_file <- file.path(map_folder, paste0(details_tab[i, 'selfmap_prefix'], ".rnaseq_gea.txt.gz"))
    # Check for results and mapping file
    if (!file.exists(hmm_file) || !file.exists(map_file)) {
      message("### Skipping ", sample_name, " for ", OG, " (file missing)")
      next  # But result_list still has 0 as default
    }
    # Load and process mapping file
    map <- read.table(map_file, sep="\t", header=TRUE,
                      colClasses=c("character", rep("NULL",5), "numeric", "numeric", 
                                   rep("NULL", 3), "numeric", rep("NULL", 3)))
    colnames(map) <- c("transcriptId", "length", "count", "countA")
    map$rpk <- map$countA / (map$length / 1000)
    map$tpm <- (map$rpk / sum(map$rpk)) * 1e6
    
    # Read and parse hmmsearch domtblout
    hmm_lines <- readLines(gzfile(hmm_file))
    hmm_data <- hmm_lines[!grepl("^#", hmm_lines)]
    if (length(hmm_data) == 0) {
      message("### Empty HMM file for ", sample_name, " in ", OG)
      next
    }
    hmm_split <- strsplit(hmm_data, "\\s+")
    hmm_df <- data.frame(
      transcriptId = sapply(hmm_split, `[`, 1),
      evalue = as.numeric(sapply(hmm_split, `[`, 5)),
      bitscore = as.numeric(sapply(hmm_split, `[`, 6)),
      stringsAsFactors = FALSE
    )
    
    # Filter
    hmm_df <- subset(hmm_df, evalue <= eval_cutoff & bitscore >= bitscore_thresh)
    if (nrow(hmm_df) == 0) {
      message("### All hits filtered out for ", sample_name, " in ", OG)
      next
    }
    
    # Join and get TPM sum
    hmm_map <- map[map$transcriptId %in% hmm_df$transcriptId, ]
    tpm_sum <- sum(hmm_map$tpm)
    # Store result (overwrite 0)
    result_list[[sample_name]] <- tpm_sum
  }
  
  # Convert result list to data frame and write to CSV
  og_result <- data.frame(sample=names(result_list), TPM=unlist(result_list))
  write.csv(og_result, file=paste0(out_folder, "/", OG, "_TPM_results.csv"), row.names=FALSE)
}
