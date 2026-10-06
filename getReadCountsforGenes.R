# Adding read counts from exon to gene level
library(dplyr)

masterFile = "/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/gene_exons/SNPs_present_in_exons.txt"
master <- read.table(masterFile, sep = "\t", header = T, as.is = T) 

for (i in 3:8){
  countsFile = paste0("/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/readCountsAtSNP/SCR00",i,"_readCountsAtSNP.txt")
  counts <- read.table(countsFile, sep = "\t", header = T, as.is = T) 
  colnames(counts) = gsub("SNP.ID", "snpID", colnames(counts))
  
  gene_counts_merge = merge(x=master, y=counts, by="snpID")
  
  gene_counts_summary <- gene_counts_merge %>%
    group_by(ENTREZID, SYMBOL) %>% 
    summarise(B6_reads = sum(B6_reads, na.rm = TRUE),
              Spret_reads = sum(Spret_reads, na.rm = TRUE))
  
  write.table(gene_counts_summary, paste0("/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/readCountsAtGene/SCR00",i,"_readCountsAtGene.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

}
