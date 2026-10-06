# This script counts reads assigned to genome1 and genome2 at SNP positions
# library(parallel)

setwd("/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/")
snp_bed <- read.table("/home/Apps/SNPsplit/refFiles/mm39/all_SNPs_SPRET_EiJ_GRCm39_header.txt", 
                      sep = "\t", header = TRUE, as.is = TRUE)
# snp_bed <- snp_bed[-c(5)]
# colnames(snp_bed) <- c("SNP-ID",	"Chromosome",	"Position",	"Strand",	"Ref/SNP")

for (i in 2:8){
# i = 1
  sample = paste0("SCR00",i)
  projectDir <- "/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/"
  bamDir <- "/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/bams/"
  B6_bam <- paste0(bamDir,sample,".sorted.genome1.bam.dupMark.bam")
  Spret_bam <- paste0(bamDir,sample,".sorted.genome2.bam.dupMark.bam")
  
  # Column names for counts
  B6colname <- "B6_reads"
  Spretcolname <- "Spret_reads"
  
  snp_bed[, B6colname] <- NA
  snp_bed[, Spretcolname] <- NA
  
  for (j in 1:nrow(snp_bed)){
    # for (j in 1:10){  # debug
    chr = snp_bed$Chromosome[j]
    start = snp_bed$Position[j]
    end = snp_bed$Position[j]
    
    print(paste(chr,":", start,"-",end))
    
    # Calculate progress percentage
    progress <- (j / nrow(snp_bed)) * 100
    timestamp <- Sys.time()
    print(sprintf("Progress: %.2f%% completed at %s", progress, timestamp))
    
    # B6 counts
    B6cmd = paste("samtools view -c -F 4", B6_bam, paste0(chr,":",start,"-", end))
    B6Output <- capture.output(system(B6cmd, intern = TRUE))
    numericValue <- as.numeric(gsub("^\\[1\\] \"(\\d+)\"$", "\\1", B6Output))
    snp_bed[[B6colname]][j] = numericValue
    
    # Spret counts
    Spretcmd = paste("samtools view -c -F 4", Spret_bam, paste0(chr,":",start,"-", end))
    spretOutput <- capture.output(system(Spretcmd, intern = TRUE))
    numericValue <- as.numeric(gsub("^\\[1\\] \"(\\d+)\"$", "\\1", spretOutput))
    snp_bed[[Spretcolname]][j] = numericValue
    
  }
  write.table(snp_bed, file = paste0(sample,"_readCountsAtSNP.txt"), sep = "\t", quote = F, row.names = F)
  
}
