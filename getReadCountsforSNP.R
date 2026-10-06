library(parallel)
setwd("/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/")

# set sample id
i = 6
sample <- paste0("SCR00", i)

# Make a log file for SANITY
log_file <- paste0(sample,"_progress_log.txt")
write(paste0("Starting job... ", timestamp <- Sys.time(), "\n"), file = log_file)

# set file paths
projectDir <- "/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/"
bamDir <- "/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/bams/"
B6_bam <- paste0(bamDir, sample, ".sorted.genome1.bam.dupMark.bam")
Spret_bam <- paste0(bamDir, sample, ".sorted.genome2.bam.dupMark.bam")

# Read in SNP file
write("Reading in SNP file...", file = log_file, append = TRUE)
snp.start = Sys.time()
snp_bed <- read.table("/home/Apps/SNPsplit/refFiles/mm39/all_SNPs_SPRET_EiJ_GRCm39_header.txt",  
                      sep = "\t", header = TRUE, as.is = TRUE)
write("Done reading SNP file.\n", file = log_file, append = TRUE)
snp.end = Sys.time()
snpTime = format(round(difftime(snp.end,snp.start,units = "mins"),2))
write(paste0("Time taken to read SNP file =", snpTime, "\n"), file = log_file, append = TRUE)

# Initialize cluster
write("Making clusters.\n", file = log_file, append = TRUE)
numCores <- 30  # Use all cores except one
cl <- makeCluster(numCores)
clusterEvalQ(cl, library(parallel))  # Load required libraries on workers

# Column names for counts
B6colname <- "B6_reads"
Spretcolname <- "Spret_reads"

snp_bed[, B6colname] <- NA
snp_bed[, Spretcolname] <- NA

# Export required variables to workers
clusterExport(cl, varlist = c("snp_bed", "B6_bam", "Spret_bam", "log_file"))

# Parallel loop with tryCatch
tryCatch({
  write("Starting parallel loop...\n", file = log_file, append = TRUE)
  snp_bed[, c(B6colname, Spretcolname)] <- t(parSapply(cl, 1:nrow(snp_bed), function(j) {
    chr <- snp_bed$Chromosome[j]
    start <- snp_bed$Position[j]
    end <- snp_bed$Position[j]
    
    print(paste(chr, ":", start, "-", end))
    
    # Calculate progress percentage
    progress <- (j / nrow(snp_bed)) * 100
    start.timestamp <- Sys.time()
    progress_message <- sprintf("Progress: %.2f%% completed at %s", progress, start.timestamp)
    write(progress_message, file = log_file, append = TRUE)
    
    # Commands for counting reads
    B6cmd <- paste("samtools view -c -F 4", B6_bam, paste0(chr, ":", start, "-", end))
    Spretcmd <- paste("samtools view -c -F 4", Spret_bam, paste0(chr, ":", start, "-", end))
    
    # Run commands and parse results
    B6_reads <- as.numeric(system(B6cmd, intern = TRUE))
    Spret_reads <- as.numeric(system(Spretcmd, intern = TRUE))
    
    return(c(B6_reads, Spret_reads))
  }))
  
  end.timestamp <- Sys.time()
  write(paste0("Ending the parallel job: ", end.timestamp, "\n"), file = log_file, append = TRUE)
  # timeTaken = format(round(difftime(end.timestamp,start.timestamp, units = "hours"),2))
  write(paste0("Final timestamp:", end.timestamp, "\n"), file = log_file, append = TRUE)
  
  # Write results for the current sample
  write("Writing read counts file...\n", file = log_file, append = TRUE)
  write.table(snp_bed, file = paste0("/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/readCountsAtSNP/", sample, "_readCountsAtSNP.txt"), sep = "\t", quote = F, row.names = F)
  
}, error = function(e) {
  # This block catches errors and prints the message to the log file
  cat("Error:", e$message, "\n")
  write(paste("Error:", e$message), file = log_file, append = TRUE)
  
}, finally = {
  # Always stop the cluster, even if an error occurs
  stopCluster(cl)
  write("All done! Stopped the clusters!", file = log_file, append = TRUE)
})

# Log the overall time taken
write(paste0("Total Time taken:", format(round(difftime(Sys.time(),snp.start, units = "days"),2)), "\n"), file = log_file, append = TRUE)
