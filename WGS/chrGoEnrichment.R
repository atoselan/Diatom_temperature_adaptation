# Import libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("topGO")
library(topGO)
library(ggplot2)


# Load GO mapping for entire genome (our background)
go_data <- read.table("C://Users/andre/OneDrive/Documents/Work/EvoExp/rnaSeq/TP_data/Thaps3_go_for_trinity_de.txt",
                     sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(go_data) <- c("gene", "GO")
   
geneID2GO <- strsplit(go_data$GO, split=",")
names(geneID2GO) <- go_data$gene
   
# Load target genes (genes to test)
target_genes <- read.table("C://Users/andre/OneDrive/Documents/Work/EvoExp/rnaSeq/TP_data/TP_DE_3/gene_group_expression/chr_14/chr14_ids.txt",
                          header=FALSE, stringsAsFactors=FALSE)$V1
target_genes <- gsub("^\\^", "", target_genes)
   
# Create named vector of all genes (0 = background, 1 = target)
allGenes <- factor(as.integer(names(geneID2GO) %in% target_genes))
names(allGenes) <- names(geneID2GO)


# Extract genes annotated to each significant GO term
get_genes_for_terms <- function(GOdata, sigTerms){
  gene_list <- lapply(sigTerms$GO.ID, function(go_id){
    # all genes annotated to this GO term
    all_genes <- genesInTerm(GOdata, go_id)[[1]]
    # subset to target genes only
    target_hits <- all_genes[all_genes %in% names(allGenes[allGenes == 1])]
    data.frame(GO.ID = go_id,
               Term = sigTerms$Term[sigTerms$GO.ID == go_id],
               TargetGenes = paste(target_hits, collapse=","),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, gene_list)
}


# Function to run topGO for one ontology
run_topGO <- function(ont){
  
  GOdata <- new("topGOdata",
                ontology = ont,
                allGenes = allGenes,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO)
     
  resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
     
  allRes <- GenTable(GOdata, classicFisher=resultFisher,
                     orderBy="classicFisher", ranksOf="classicFisher", topNodes=100)
     
  # Adjust p-value for multiple testing (BH)
  allRes$FDR <- p.adjust(allRes$classicFisher, method="fdr")
     
  # Keep significant terms
  sigGO <- allRes[allRes$FDR <= 0.05, ]
  
  # Get genes from significant GO terms
  get_genes_for_terms(GOdata, sigGO)
  
  if(nrow(sigGO)==0) return(NULL)
    sigGO <- sigGO %>%
      mutate(logFDR = -log10(FDR),
             Ontology = ont,
             Term = factor(Term, levels=rev(Term))) # for plotting
     
  return(sigGO)
}


# Run GO enrichment for BP, MF, CC
sigBP <- run_topGO("BP")
sigMF <- run_topGO("MF")
sigCC <- run_topGO("CC")

# Combine results
sigAll <- bind_rows(sigBP, sigMF, sigCC)
   
# Set colors for bars
bar_colors <- c(MF="#33a02c", BP="#1f78b4", CC="#e31a1c")

# Make faceted barplot
p <- ggplot(sigAll, aes(x=Term, y=logFDR, fill=Ontology)) +
  geom_col(color = "black", size = 0.3) +
  coord_flip() +
  facet_grid(
    rows = vars(Ontology), 
    scales = "free_y", 
    space = "free", 
    switch = "x"
  ) +
  scale_fill_manual(values = bar_colors) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background.y = element_rect(fill="#f0f0f0", color="black", size=0.8),
    strip.text.y.right = element_text(
    angle = 270,        # rotate text 90 degrees
    face = "bold",     # bold
    color = "black",
    size = 10,
    hjust = 0.5,       # adjust horizontal position
    vjust = 0.5,       # adjust vertical position
    margin = margin(r = 5, l = 5,
                    t = 5, b = 5)
    ),
    panel.border = element_rect(color="black", fill=NA, size=0.8),
    panel.spacing = unit(0.5, "lines"),
    axis.text.y = element_text(size=10, face="bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color="grey80"),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "", y = "-log10(FDR)", fill = "Ontology")


# Save plot to png
ggsave(
  filename = "C://Users/andre/OneDrive - University of East Anglia/EvoExp_notes/FC_patched_work/Allele_balance/chr12_chr14_GOenrich/chr_12__enrichedGO_FDR_05.png",
  plot = p,
  bg="white",
  width = 8,
  height = 8,
  dpi = 300
)
