# Import libraries
library("vcfR")

# Set working directory
setwd("C://Users/andre/OneDrive/Documents/Work/EvoExp/ShinyNewSnps")

# Avoid scientific notation for genome position
options(scipen=999)

# Read in variant file, contains both snps and indels
vcf_file <- read.vcfR("chrSpecificVcf/evoExp.chr_14.bcfCall.vars.vcf.gz",
                      verbose=FALSE)

# Remove samples predicted to be outliers, replicates
# 1 and 4 of T4_32C
# These are 'columns' 39 and 42 of the vcf
tmp <- vcf_file[,c(-39,-42)]
vcf_file <- tmp
rm(tmp)

# Filter out indels
snpOnly <- extract.indels(vcf_file)

# Remove non bi-allelic loci
bi <- is.biallelic(snpOnly)
biSnpOnly <- snpOnly[bi,]

# Delete intermediate files from memory
rm(bi, vcf_file, snpOnly)

# Get allele depths
ad <- extract.gt(biSnpOnly, element = 'AD')
# Get allele specific depths
allele1 <- masplit(ad, record = 1)
allele2 <- masplit(ad, record = 2)
# Get allele %
ad1 <- allele1 / (allele1 + allele2)
ad2 <- allele2 / (allele1 + allele2)
# Get allele frequencies
freq1 <- ad1/(ad1+ad2)
freq2 <- ad2/(ad1+ad2)

# Find density peaks
myPeaks1 <- freq_peak(freq1, getPOS(biSnpOnly), winsize = winsize, bin_width = bin_width)
myPeaks2 <- freq_peak(freq2, getPOS(biSnpOnly), winsize = winsize, bin_width = bin_width, lhs = FALSE)

# Filter out low count peaks
is.na(myPeaks1$peaks[myPeaks1$counts < 20]) <- TRUE
is.na(myPeaks2$peaks[myPeaks2$counts < 20]) <- TRUE

# Set plotting parameters
winsize <- 1e5
bin_width <- 0.02
sampleNames <- c("T0_22C", "T2_32C")

# Plot allele balance data
pdf("alleleBalancePlots/chr_14_T0-T2_alleleBalance.pdf", height=6, width=10)
  # Matrix layout for 2 samples
  layout(matrix(c(1,3,2,4), nrow=2, ncol=2), 
         widths=c(4,1,4,1), 
         heights=c(3,3))
  j <- 1
  # Loop through samples of interest 
  for (i in c(2,24)){    
    par(mar=c(3,4,3,0))
    # Get sample name
    mySample <- colnames(freq1)[i]
    plot(getPOS(biSnpOnly), freq1[,mySample], ylim=c(0,1), type="n", yaxt='n', 
         main = paste(sampleNames[j], ": chr_14", sep=''), xlab = "POS", ylab="Allele balance")
    j <- j+1
    axis(side=2, at=c(0,0.25,0.333,0.5,0.666,0.75,1), 
         labels=c(0,'1/4','1/3','1/2','2/3','3/4',1), las=1)
    abline(h=c(0.25,0.333,0.5,0.666,0.75), col=8)
    points(getPOS(biSnpOnly), freq1[,mySample], pch=20, col="#A6CEE344")
    points(getPOS(biSnpOnly), freq2[,mySample], pch=20, col="#1F78B444")
    segments(x0=myPeaks1$wins[,'START_pos'], y0=myPeaks1$peaks[,mySample],
             x1=myPeaks1$wins[,'END_pos'], lwd=2)
    segments(x0=myPeaks1$wins[,'START_pos'], y0=myPeaks2$peaks[,mySample],
             x1=myPeaks1$wins[,'END_pos'], lwd=2)
    # Prep histogram data
    bp1 <- hist(freq1[,mySample], breaks=seq(0,1,by=bin_width), plot=FALSE)
    bp2 <- hist(freq2[,mySample], breaks=seq(0,1,by=bin_width), plot=FALSE)
    # Plot histogram
    par(mar=c(3,1,3,2))
    barplot(height=bp1$counts, width=0.02,  space=0, horiz=T, add=FALSE, col="#A6CEE3")
    barplot(height=bp2$counts, width=0.02,  space=0, horiz=T, add=TRUE, col="#1F78B4")
  }
dev.off()
