# -------------------------------
# Import libraries
# -------------------------------
# BiocManager::install("SNPRelate")
#BiocManager::install("snpStats")
library(SNPRelate)
library(vcfR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(tidyr)
library(pheatmap)
library(inflection)
library(topGO)
library(Polychrome)


# -------------------------------
# Basic setup
# -------------------------------
# Set working directory
setwd("C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/TP_chromosomal_stability_genes/SnpEff_approach")

# Set file paths
vcf_file <- "TP_vars_biallelic_snps.vcf.gz"
gds_file <- "filtered_snps.gds"
temp_vcf <- "filtered_temp.vcf"

# Abbreviated sample names - change old names to these for plotting
new_sample_names <- c("T0_22C_1", "T0_22C_2", "T0_22C_3",
                      "T3_32C_1",
                      "T1_22C_1", "T1_22C_2", "T1_22C_3", "T1_22C_4", "T1_22C_5",
                      "T1_32C_2", "T1_32C_3", "T1_32C_4", "T1_32C_5",
                      "T1_9C_1", "T1_9C_2", "T1_9C_3", "T1_9C_4", "T1_9C_5",
                      "T2_22C_2", "T2_22C_3", "T2_22C_4", "T2_22C_5","T2_32C_1", 
                      "T2_32C_2", "T2_32C_3", "T2_32C_4", "T2_32C_5",
                      "T2_9C_1", "T2_9C_2", "T2_9C_3", "T2_9C_4", "T2_9C_5",
                      "T3_9C_1", "T3_9C_2", "T3_9C_3", "T3_9C_4", "T3_9C_5", 
                      "T4_32C_1", "T4_32C_2", "T4_32C_3", "T4_32C_4", "T4_32C_5")


# -------------------------------
# Read and rename VCF
# -------------------------------
vcf <- read.vcfR(vcf_file, verbose=FALSE)
#table(extract.gt(vcf, element="GT"))
# 0/0     0/1     0/2     0/3     1/1     1/2     1/3     2/.     2/0     2/1     2/2     2/3     3/1     3/3 
# 190592 7900396    2773      17   58499     664       1    1661     179     344    1462       6       2       1 

# Extract chromosome names from header
chrom_order <- sub(".*ID=([^,>]+).*", "\\1", grep("^##contig", vcf@meta, value = TRUE))

# Rename samples 
colnames(vcf@gt)[-1] <- new_sample_names


# -------------------------------
# Filter VCF by depth, genotype quality, and missingness
# -------------------------------
# Extract genotype, depth, and GQ
gt <- extract.gt(vcf, element="GT", as.numeric=FALSE)
dp <- extract.gt(vcf, element="DP", as.numeric=TRUE)
gq <- extract.gt(vcf, element="GQ", as.numeric=TRUE)

gt_masked <- gt
gt_masked[dp < 10] <- NA  # Mask low depth genotypes
gt_masked[gq < 30] <- NA  # Mask low genotype quality genotypes
gt_masked[gt_masked %in% c("./.", ".")] <- NA

# Check missingness per sample
# Note - this is the mean of a boolean
missing_per_sample <- colMeans(is.na(gt_masked))
#print(missing_per_sample)
barplot(sort(missing_per_sample),
        las = 2, cex.names = 0.6,
        ylab = "Fraction missing genotypes",
        ylim=c(0,0.2),
        main = "Per-sample missingness")


# Remove SNPs with >10% missing genotypes
missing_per_snp <- rowMeans(is.na(gt_masked))
keep_snps <- which(missing_per_snp <= 0.1)
vcf_filtered <- vcf[keep_snps, ]
table(extract.gt(vcf_filtered, element="GT"))

# Write filtered VCF to temporary file
write.vcf(vcf_filtered, file=temp_vcf)

# Tidy up - remove large files no longer needed
rm(vcf, gt, gt_masked, dp, gq, keep_snps, missing_per_snp)
gc()


# -------------------------------
# Run PCA
# -------------------------------
# Convert filtered vcf to GDS format for SNPRelate
snpgdsVCF2GDS(temp_vcf, gds_file, method="biallelic.only")
genofile <- snpgdsOpen(gds_file)

# Get snp specific info - useful for indexing later
# by chromosome and position rather than snp_id from
# the PCA loadings
gdssnp_info <- snpgdsSNPList(genofile)
gdssnp_info$key <- paste0(gdssnp_info$chromosome, ":", gdssnp_info$position)

# Define samples to include in PCA
sample_list <- c("T0_22C_1", "T0_22C_2", "T0_22C_3",
                 "T1_22C_1", "T1_22C_2", "T1_22C_3", "T1_22C_4", "T1_22C_5",
                 "T1_32C_2", "T1_32C_3", "T1_32C_4", "T1_32C_5",
                 "T1_9C_1", "T1_9C_2", "T1_9C_3", "T1_9C_4", "T1_9C_5",
                 "T2_22C_2", "T2_22C_3", "T2_22C_4", "T2_22C_5","T2_32C_1", 
                 "T2_32C_2", "T2_32C_3", "T2_32C_4", "T2_32C_5",
                 "T3_32C_1",
                 "T2_9C_1", "T2_9C_2", "T2_9C_3", "T2_9C_4", "T2_9C_5",
                 "T3_9C_1", "T3_9C_2", "T3_9C_3", "T3_9C_4", "T3_9C_5", 
                 "T4_32C_2", "T4_32C_3", "T4_32C_5")

# Run PCA
pca <- snpgdsPCA(genofile, autosome.only=FALSE, sample.id = sample_list)
pc_percent <- pca$varprop * 100  # % variance explained

# PCA dataframe for plotting 1st 2 PCs
pca_df <- data.frame(
  sample = pca$sample.id,
  PC1 = pca$eigenvect[,1],
  PC2 = pca$eigenvect[,2]
)

# Color by temperature group
pca_df$group <- ifelse(grepl("32C", pca_df$sample), "32C",
                       ifelse(grepl("22C", pca_df$sample), "22C",
                              ifelse(grepl("9C", pca_df$sample), "9C", "other")))
group_colors <- c("32C" = "red", "22C" = "green", "9C" = "blue", "other" = "gray")

# Plot PCA
pca_p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = group), color = "black", shape = 21, size = 3, stroke = 0.7) +
  stat_ellipse(aes(fill = group), geom = "polygon", type = "norm", level = 0.95,
               alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(values = group_colors) +
  geom_text_repel(data = subset(pca_df, group == "32C"), 
                  aes(label = sample), size = 2, max.overlaps = Inf) +
  labs(x=paste0("PC1 (", round(pc_percent[1],1), "%)"),
       y=paste0("PC2 (", round(pc_percent[2],1), "%)"),
       title="PCA of filtered SNPs") +
  theme_bw() +
  theme(panel.grid = element_blank())

# Save as high-quality PNG
ggsave("PCA_snpBased_plot.png", plot = pca_p,
       width = 4, height = 4, units = "in", dpi = 600)


# -------------------------------
# Extract SNP loadings
# -------------------------------
# Principal component loadings are the coefficients (weights)
# for each SNP, showing how much it contributes to a given PC. 
# A PC is a linear combination of the original SNPs, weighted 
# by these loadings.
snp_loadings <- snpgdsPCASNPLoading(pca, genofile)

# Close GDS - don't use this anymore
snpgdsClose(genofile)
rm(genofile)
gc()

# Get SNP position data from GDS
pca_snp_info <- gdssnp_info[snp_loadings$snp.id, c("snp.id","chromosome","position","key")]
#head(pca_snp_info)

# Create keys for the filtered VCF
vcf_keys <- paste0(vcf_filtered@fix[,"CHROM"], ":", vcf_filtered@fix[,"POS"])
vcf_keys <- gsub("chr", "", vcf_keys)                    

# Make a vector of INFO
vcf_info <- vcf_filtered@fix[,"INFO"]

# Create a named vector for easy lookup
names(vcf_info) <- vcf_keys

# Now for your PCA SNPs:
pca_snp_info$INFO <- vcf_info[pca_snp_info$key]

# Combine loadings with SNP info
loading_df <- data.frame(
  SNP = pca_snp_info$snp.id,
  CHR = pca_snp_info$chromosome,
  POS = pca_snp_info$position,
  INFO = pca_snp_info$INFO,
  PC1_loading = snp_loadings$snploading[1, ]  # first row = PC1
)

rm(pca_snp_info, vcf_info)
gc()


# -------------------------------
# Filter PCA loadings by inflection point
# -------------------------------
get_top_snps <- function(loading_df, pc_col = "PC1_loading", sign = c("pos","neg")) {
  
  sign <- match.arg(sign)
  
  if(sign == "pos") {
    snps <- loading_df[loading_df[[pc_col]] > 0, ]
    snps <- snps[order(-snps[[pc_col]]), ]
  } else {
    snps <- loading_df[loading_df[[pc_col]] < 0, ]
    snps <- snps[order(-snps[[pc_col]]), ] # ascending for negative
  }
  
  x <- seq_len(nrow(snps))
  y <- snps[[pc_col]]
  ip <- findiplist(x, y, index = 1)
  
  # Choose cutoff
  cutoff_idx <- ip["EDE","j2"]
  # Test with no cutoff
  #cutoff_idx <- NROW(snps)
  top_snps <- head(snps, cutoff_idx)
  
  # Optional: plot
  png(paste0("PC1_", sign, "_loadings_elbow_plot.png"),
      height=1200, width=1500, res=300)
  plot(x, y, type="l",
       main = paste0("Elbow points for ", sign, " SNP loadings"),
       xlab="Index score", ylab=paste0(pc_col, " score"))
  abline(v = ip["EDE","j1"], col = "blue", lty=2, lwd=2)
  abline(v = ip["EDE","j2"], col = "red", lty=2, lwd=2)
  legend("topright",
         legend = c("j1","j2"),
         col = c("blue","red"),
         lty=2, lwd=2, cex=0.8)
  dev.off()
  
  return(top_snps)
}

top_positive_snps <- get_top_snps(loading_df, sign="pos")
top_negative_snps <- get_top_snps(loading_df, sign="neg")


# -------------------------------
# Annotate top SNPs from snpEff INFO
# -------------------------------
# Function to annotate top SNPs from snpEff INFO field
annotate_snps <- function(top_snps) {
  
  # Internal helper to parse a single INFO string
  parse_ann <- function(info_str) {
    ann_field <- str_match(info_str, "ANN=([^;]+)")[,2]
    if (is.na(ann_field)) return(NA)
    # Take first annotation if multiple
    ann_main <- strsplit(ann_field, ",")[[1]][1]
    # Split by | according to snpEff format
    ann_parts <- unlist(strsplit(ann_main, "\\|"))
    # Build a small tibble
    tibble(
      allele     = ann_parts[1],
      effect     = ann_parts[2],
      impact     = ann_parts[3],
      gene_id    = ann_parts[6],
      nt_change  = ann_parts[12],
      aa_change  = ann_parts[13]
    )
  }
  
  # Apply parsing to all SNPs
  parsed_ann <- lapply(top_snps$INFO, parse_ann)
  parsed_ann <- bind_rows(parsed_ann)
  
  # Combine with SNP metadata
  snps_annotated <- cbind(top_snps, parsed_ann)
  
  # Create a key
  snps_annotated$key <- paste0(snps_annotated$CHR, ":", snps_annotated$POS)
  
  # Drop the raw INFO column
  snps_annotated$INFO <- NULL
  
  return(snps_annotated)
}

top_positive_snps_annotated <- annotate_snps(top_positive_snps)
top_negative_snps_annotated <- annotate_snps(top_negative_snps)

# Clean up original top SNP objects
rm(top_positive_snps, top_negative_snps)
gc()


# -------------------------------
# Add genotypes for top SNPs
# -------------------------------
# Extract genotype matrix from filtered VCF
gt_matrix <- extract.gt(vcf_filtered, element="GT", as.numeric=FALSE)

# Match positions based on genome positional key
snp_idx_p <- match(top_positive_snps_annotated$key, vcf_keys)
snp_idx_n <- match(top_negative_snps_annotated$key, vcf_keys)

# Check for failed matches
sum(is.na(snp_idx_p))
sum(is.na(snp_idx_n))

gt_for_top_pos <- gt_matrix[snp_idx_p, , drop = FALSE]
gt_for_top_neg <- gt_matrix[snp_idx_n, , drop = FALSE]
# combine with top positive SNP table
top_with_gt_pos <- cbind(top_positive_snps_annotated, as.data.frame(gt_for_top_pos))
top_with_gt_neg <- cbind(top_negative_snps_annotated, as.data.frame(gt_for_top_neg))

# Clean up
rm(gt_matrix, snp_idx_p, snp_idx_n, gt_for_top_pos, gt_for_top_neg)
gc()


# -------------------------------
# Haplotype block calculation
# -------------------------------
haplotype_blocks <- function(snp_table, dist_cutoff = 10000, r2_cutoff = 0.5, min_samples_frac = 0.9) {
  
  # Extract genotype matrix (assumes sample columns start with "T")
  geno_mat <- as.matrix(snp_table[, grep("^T", colnames(snp_table))])
  
  # Convert to numeric allele counts
  geno_num <- matrix(NA_integer_, nrow = nrow(geno_mat), ncol = ncol(geno_mat),
                     dimnames = dimnames(geno_mat))
  geno_num[geno_mat == "0/0"] <- 0L
  geno_num[geno_mat %in% c("0/1", "1/0")] <- 1L
  geno_num[geno_mat == "1/1"] <- 2L
  geno_num[is.na(geno_mat) | geno_mat == "./."] <- NA_integer_
  
  # Assign blocks on a single chromosome
  assign_blocks_chr_numeric <- function(chr, pos, geno_chr_num) {
    n <- length(pos)
    if (n == 1) return(data.frame(CHR = chr, POS = pos, block_id = 1L))
    
    block_id <- rep(NA_integer_, n)
    block <- 1L
    i <- 1L
    min_samples <- ceiling(min_samples_frac * ncol(geno_chr_num))
    
    while (i <= n) {
      j <- i
      block_snps <- i  # indices of SNPs currently in block
      
      repeat {
        # candidate SNPs within dist_cutoff of the furthest SNP in the block
        block_end_pos <- max(pos[block_snps])
        next_snp <- j + 1
        if (next_snp > n) break
        if ((pos[next_snp] - block_end_pos) > dist_cutoff) {
          j <- next_snp - 1
          break
        }
        
        # r² between candidate SNP and all SNPs in the current block
        ld_vals <- sapply(block_snps, function(k) {
          pair_idx <- which(!is.na(geno_chr_num[k, ]) & !is.na(geno_chr_num[next_snp, ]))
          if (length(pair_idx) < min_samples) return(NA_real_)
          x <- geno_chr_num[k, pair_idx]
          y <- geno_chr_num[next_snp, pair_idx]
          if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
          cor(x, y)^2
        })
        
        ld_ok <- all(!is.na(ld_vals) & ld_vals > r2_cutoff)
        if (!ld_ok) break
        
        # extend block
        j <- next_snp
        block_snps <- c(block_snps, next_snp)
      }
      
      block_id[i:j] <- block
      block <- block + 1L
      i <- j + 1L
    }
    data.frame(CHR = chr, POS = pos, block_id = block_id)
  }
  
  # Split SNPs by chromosome
  chr_split <- split(seq_len(nrow(snp_table)), snp_table$CHR)
  
  # Apply per chromosome
  all_blocks <- do.call(rbind, lapply(chr_split, function(idx) {
    chr <- snp_table$CHR[idx][1]
    pos <- snp_table$POS[idx]
    geno_chr_num <- geno_num[idx, , drop = FALSE]
    assign_blocks_chr_numeric(chr, pos, geno_chr_num)
  }))
  
  # Merge block info back into snp_table
  snp_table$block_id <- all_blocks$block_id
  snp_table <- snp_table %>%
    group_by(CHR, block_id) %>%
    mutate(snps_in_block = n(),
           block_start = min(POS),
           block_end   = max(POS),
           block_size = (block_end - block_start)) %>%
    ungroup()
  
  # Split table for convenience
  list(
    all = snp_table,
    filtered = snp_table %>% filter(snps_in_block <= 2 & block_size < 100),   # isolated SNPs
    blocks   = snp_table %>% filter(snps_in_block > 2 & block_size > 100)     # LD blocks
  )
}

# Run for positive SNPs
haplo_pos <- haplotype_blocks(top_with_gt_pos)
# Run for negative SNPs
haplo_neg <- haplotype_blocks(top_with_gt_neg)


# -------------------------------
# Plot genotype heatmaps
# -------------------------------
plot_snp_heatmap <- function(dataset, chrom_order, sample_order, 
                             block_filter = c("all", "filtered", "blocks"),
                             out_file = NULL) {
  block_filter <- match.arg(block_filter)
  
  # Subset depending on block filter
  if (block_filter == "filtered") {
    df <- as.data.frame(dataset[dataset$snps_in_block <= 2, ])
  } else if (block_filter == "blocks") {
    df <- as.data.frame(dataset[dataset$snps_in_block > 2, ])
  } else {
    df <- as.data.frame(dataset)
  }
  
  if (nrow(df) == 0) {
    stop("No SNPs available for the chosen block_filter")
  }
  
  # Chromosome order
  #chrom_order <- gsub("chr", "", chrom_lengths$CHR)
  df$CHR <- factor(df$CHR, levels = chrom_order)
  
  # Extract genotype columns
  gt <- df %>% dplyr::select(starts_with("T"))
  
  # Convert genotypes to numeric codes
  geno_mat <- apply(as.matrix(gt), c(1,2), function(x) {
    if (is.na(x)) return(0)   # NA → 0
    if (x == "0/0") return(1) # homozygous reference
    if (x == "0/1") return(2) # heterozygous
    if (x == "1/1") return(3) # homozygous alternate
    return(0)
  })
  
  # Row/col names
  rownames(geno_mat) <- df$key
  colnames(geno_mat) <- colnames(gt)
  
  # Sort SNPs by chromosome + position
  # Sort SNPs by chromosome + position using chrom_order
  df_sorted <- df %>%
    mutate(CHR = factor(CHR, levels = chrom_order)) %>%
    arrange(CHR, POS)
  
  chrom_present <- unique(df_sorted$CHR)
  chrom_order_filtered <- chrom_order[chrom_order %in% chrom_present]
  
  geno_mat <- geno_mat[df_sorted$key, , drop = FALSE]
  
  # Chromosome annotation for rows
  row_anno <- df_sorted %>%
    dplyr::select(CHR) %>%
    mutate(CHR = factor(CHR, levels = chrom_order_filtered))
  rownames(row_anno) <- df_sorted$key
  
  # Chromosome colors
  chrom_cols <- setNames(
    createPalette(length(chrom_order_filtered), 
                  c("#ff0000", "#00ff00", "#0000ff")),
    chrom_order_filtered
  )
  ann_colors <- list(CHR = chrom_cols)
  
  # Genotype colors
  my_colors <- c("white", "#4575b4", "#fee090", "#a50026")
  
  # Heatmap
  p <- pheatmap(
    geno_mat[, sample_order, drop = FALSE],
    color = my_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_anno,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    legend_breaks = 0:3,
    legend_labels = c("NA", "0/0", "0/1", "1/1")
  )
  
  # Save if requested
  if (!is.null(out_file)) {
    ggsave(out_file, plot = p, width = 8, height = 6, dpi = 300)
  }
  
  return(p)
}

# Sample order for plotting
sample_order <- c("T0_22C_1", "T0_22C_2", "T0_22C_3",
                  "T1_22C_1", "T1_22C_2", "T1_22C_3", "T1_22C_4", "T1_22C_5",
                  "T2_22C_2", "T2_22C_3", "T2_22C_4", "T2_22C_5", 
                  "T1_9C_1", "T1_9C_2", "T1_9C_3", "T1_9C_4", "T1_9C_5",
                  "T2_9C_1", "T2_9C_2", "T2_9C_3", "T2_9C_4", "T2_9C_5",
                  "T3_9C_1", "T3_9C_2", "T3_9C_3", "T3_9C_4", "T3_9C_5", 
                  "T1_32C_2", "T1_32C_3", "T1_32C_4", "T1_32C_5",
                  "T2_32C_1", "T2_32C_2", "T2_32C_3", "T2_32C_4", "T2_32C_5",
                  "T3_32C_1",
                  "T4_32C_2", "T4_32C_3", "T4_32C_5")

# Fix chromosome naming to match chrom_order
fix_chr_names <- function(chr_vec) {
  chr_vec <- gsub("^_", "chr_", chr_vec)   # add "chr_" to those starting with "_"
  chr_vec <- gsub("^chr_(\\d+)$", "chr_\\1", chr_vec) # ensure chr_ prefix
  return(chr_vec)
}

# Update chromosome names, just leave _XX
haplo_pos$all$CHR <- fix_chr_names(haplo_pos$all$CHR)
haplo_neg$all$CHR <- fix_chr_names(haplo_neg$all$CHR)

# Positive snp - blocks >2 snps per block, filtered 1 or 2
plot_snp_heatmap(haplo_pos$all, 
                 chrom_order, sample_order, 
                 block_filter = "blocks",
                 out_file = "posSnps_LDblocks_genotype_heatmap.png")
plot_snp_heatmap(haplo_pos$all, chrom_order, sample_order, 
                 block_filter = "filtered",
                 out_file = "posSnps_LDfiltered_genotype_heatmap.png")
# Negative snps - blocks >2 snps per block, filtered 1 or 2
plot_snp_heatmap(haplo_neg$all, chrom_order, sample_order, 
                 block_filter = "blocks",
                 out_file = "negSnps_LDblocks_genotype_heatmap.png")
plot_snp_heatmap(haplo_neg$all, chrom_order, sample_order, 
                 block_filter = "filtered",
                 out_file = "negSnps_LDfiltered_genotype_heatmap.png")


# -------------------------------
# Plot allele frequencies for high impact snps
# -------------------------------
plot_high_snps <- function(snp_df, vcf, sample_order, impact_filter = "HIGH") {
  
  # Filter SNPs by impact
  snp_high <- snp_df %>%
    filter(impact == impact_filter) %>%
    mutate(CHROM_POS = paste("chr", CHR, ":", POS, sep = ""),
           SNP_label = paste0("chr", CHR, ":", POS, " (", gene_id, ")"))
  
  if (nrow(snp_high) == 0) {
    stop("No SNPs found with impact = ", impact_filter)
  }
  
  # Match to VCF rows
  vcf_pos <- paste(vcf@fix[,"CHROM"], vcf@fix[,"POS"], sep=":")
  snp_indices <- which(vcf_pos %in% snp_high$CHROM_POS)
  vcf_sub <- vcf[snp_indices, ]
  
  # Extract allele depths
  ad <- extract.gt(vcf_sub, element = "AD", as.numeric = FALSE)
  
  parse_ad_robust <- function(ad_str) {
    if (is.na(ad_str) || ad_str == "" || ad_str == ".") return(c(NA, NA))
    counts <- unlist(strsplit(ad_str, ","))
    counts_num <- sapply(counts, function(x) ifelse(x %in% c(".", ""), 0, as.numeric(x)))
    if (length(counts_num) < 2) counts_num <- c(counts_num[1], 0)
    return(counts_num[1:2]) # REF, ALT
  }
  
  n_snps <- nrow(ad)
  n_samples <- ncol(ad)
  ad_counts <- array(NA, dim = c(n_snps, n_samples, 2),
                     dimnames = list(rownames(ad), colnames(ad), c("REF","ALT")))
  
  for (i in 1:n_snps) {
    for (j in 1:n_samples) {
      ad_counts[i,j,] <- parse_ad_robust(ad[i,j])
    }
  }
  
  # Calculate allele frequencies
  alt_freq <- ad_counts[,,2] / rowSums(ad_counts, dims = 2)
  rownames(alt_freq) <- snp_high$SNP_label
  alt_freq <- alt_freq[, sample_order, drop = FALSE]
  
  # Tidy data for plotting
  df <- as.data.frame(alt_freq, check.names = FALSE)
  df$SNP <- rownames(df)
  
  long <- df %>%
    pivot_longer(cols = -SNP, names_to = "Sample", values_to = "ALT_AF") %>%
    mutate(REF_AF = 1 - ALT_AF)
  
  plot_df <- long %>%
    pivot_longer(cols = c(ALT_AF, REF_AF),
                 names_to = "Allele", values_to = "AF") %>%
    mutate(Allele = recode(Allele, ALT_AF = "ALT", REF_AF = "REF"),
           Sample = factor(Sample, levels = sample_order))
  
  # Plot
  p <- ggplot(plot_df, aes(x = Sample, y = AF, color = Allele, group = Allele)) +
    geom_hline(yintercept = c(0, 0.5, 1.0), linetype = "dotted", color = "grey40") +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE, size = 1.6) +
    facet_wrap(~ SNP, ncol = 1, scales = "fixed") +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      strip.text = element_text(size = 10),
      panel.spacing = unit(0.4, "lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Sample", y = "Allele frequency")
  
  return(p)
}

p1 <- plot_high_snps(haplo_pos$blocks, vcf_filtered, sample_order, impact_filter = "HIGH")
p1 <- plot_high_snps(haplo_pos$filtered, vcf_filtered, sample_order, impact_filter = "HIGH")
print(p1)


# -------------------------------
# GO enrichment
# -------------------------------
run_go_enrichment <- function(snp_df,
                              go_mapping_file,
                              impact_filter = c("HIGH", "MODERATE"),
                              fisher_cutoff = 0.05,
                              top_nodes = 100) {
  stopifnot(file.exists(go_mapping_file))
  
  # Filter SNPs by impact
  target_snps <- snp_df %>% filter(impact %in% impact_filter)
  target_genes <- unique(target_snps$gene_id)
  if (length(target_genes) == 0) {
    warning("No genes found with given impact filter.")
    return(NULL)
  }
  
  # Load GO mapping
  go_data <- read.table(go_mapping_file,
                        sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(go_data) <- c("gene", "GO")
  
  geneID2GO <- strsplit(go_data$GO, split = ",")
  names(geneID2GO) <- go_data$gene
  
  # Background set (all genes in mapping)
  allGenes <- factor(as.integer(names(geneID2GO) %in% target_genes))
  names(allGenes) <- names(geneID2GO)
  
  # Helper: extract genes for sig terms
  get_genes_for_terms <- function(GOdata, sigTerms) {
    gene_list <- lapply(sigTerms$GO.ID, function(go_id) {
      all_genes <- genesInTerm(GOdata, go_id)[[1]]
      target_hits <- all_genes[all_genes %in% names(allGenes[allGenes == 1])]
      data.frame(GO.ID = go_id,
                 Term = sigTerms$Term[sigTerms$GO.ID == go_id],
                 TargetGenes = paste(target_hits, collapse = ","),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, gene_list)
  }
  
  # Helper: run one ontology
  run_topGO2 <- function(ont) {
    GOdata <- new("topGOdata",
                  ontology = ont,
                  allGenes = allGenes,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)
    
    resultFisher <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
    
    allRes <- GenTable(GOdata, classicFisher = resultFisher,
                       orderBy = "classicFisher",
                       ranksOf = "classicFisher",
                       topNodes = top_nodes)
    
    # Adjust FDR
    allRes$FDR <- p.adjust(allRes$classicFisher, method = "fdr")
    allRes$classicFisher <- as.numeric(allRes$classicFisher)
    
    sigGO <- allRes[allRes$classicFisher <= fisher_cutoff, ]
    if (nrow(sigGO) == 0) return(NULL)
    
    sigGO <- sigGO %>%
      mutate(logFDR = -log10(FDR),
             logFisher = -log10(classicFisher),
             Ontology = ont,
             Term = factor(Term, levels = rev(unique(Term))))
    
    # Attach gene list for each term
    genes_df <- get_genes_for_terms(GOdata, sigGO)
    sigGO <- left_join(sigGO, genes_df, by = c("GO.ID", "Term"))
    
    sigGO
  }
  
  # Run for each ontology
  sigBP <- run_topGO2("BP")
  sigMF <- run_topGO2("MF")
  sigCC <- run_topGO2("CC")
  
  sigAll <- bind_rows(sigBP, sigMF, sigCC)
  if (is.null(sigAll) || nrow(sigAll) == 0) {
    message("No significant GO terms found.")
    return(NULL)
  }
  
  # Plot
  bar_colors <- c(MF = "#33a02c", BP = "#1f78b4", CC = "#e31a1c")
  p <- ggplot(sigAll, aes(x = Term, y = logFisher, fill = Ontology)) +
    geom_col(color = "black", size = 0.3) +
    coord_flip() +
    facet_grid(rows = vars(Ontology), scales = "free_y", space = "free", switch = "x") +
    scale_fill_manual(values = bar_colors) +
    theme_minimal(base_size = 12) +
    theme(
      strip.background.y = element_rect(fill = "#f0f0f0", color = "black", size = 0.8),
      strip.text.y.right = element_text(angle = 270, face = "bold", size = 10),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text.y = element_text(size = 10, face = "bold"),
      panel.grid.major.y = element_blank()
    ) +
    labs(x = "", y = "-log10(Fisher's p)", fill = "Ontology")
  
  return(list(sig_terms = sigAll, plot = p))
}

# Define location of GO mapping file
go_file <- "C://Users/andre/OneDrive/Documents/Work/EvoExp/rnaSeq/TP_data/Thaps3_go_for_trinity_de.txt"

# Subset: HIGH + MODERATE impact SNPs in blocks
block_pos_snps <- haplo_pos$filtered %>%
  filter(impact %in% c("HIGH", "MODERATE"))

go_results <- run_go_enrichment(block_pos_snps,
                                go_mapping_file=go_file)

# Access results table
head(go_results$sig_terms)

# Plot enrichment
print(go_results$plot)

ggsave("posSnps_LDfiltered_enrichedGO.png",
       plot = go_results$plot, 
       bg="white",
       width = 8, height = 6, dpi = 300)
