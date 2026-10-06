library(ggplot2)

for (i in 1:9) {
  
  sample = paste0("SCR00",i)
  
  # Read in the data from the txt file
  data <- read.table(paste0("/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/readCountsAtSNP/",sample,"_readCountsAtSNP.txt"), header = TRUE)
  colnames(data) = gsub(paste0(sample,"_"),"", colnames(data))
  
  # Calculate the log2 ratio of B6 reads to Spret reads
  
  data$log2_ratio <- log2(data$B6_reads / data$Spret_reads)
  
  # Create histogram using ggplot2
  histogram <- ggplot(data, aes(x = log2_ratio)) +
    geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +
    labs(title = "Histogram of log2(B6 reads/Spret reads)", x = "Log2 Ratio", y = "Frequency")
  
  # Save the plot
  ggsave(paste0(sample,"_histogram.png"), plot = histogram, width = 8, height = 6, dpi = 300)
  
}


