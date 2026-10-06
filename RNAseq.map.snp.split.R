# This script maps RNA data and run SNPsplit to get parental tagged reads

# Read in bistools; sources many important functions for this script
source("/home/Apps/bitbucket/bistools/ESB_bisTools.R") # Magic
source("/home/Apps/bitbucket/sequencing-sample-qc/mapping_barplot.R") # Magic, plotting scripts
library("stringr") # str_split_1()

projectDir = "/PROJECT/DIR/" # Replace, must have '/' at the end

# Set working directory
homeDir = paste0(projectDir, "pipeline/")
setwd(homeDir)

# Manifest file
filesDir = homeDir # Replace, if necessary
filesFile = "RNAseq.sample.manifest.txt" # Replace, if necessary
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = T, as.is = T)

# Set max resources; ESB has 128 threads and 512G, allowing for multiple of these at the default of 32T/32GB to run concurrently
# Replace, adjust according to priority
threads = 32
sortMem = "32G"

# Adapter trimming options
# Replace, choose one from: nextera, illumina
# "nextera" signifies Nextera Tn5 adapters, "illumina" signifies TRUSEQ adapters
seq_platform = "nextera"

# Picard call
picardCmd = "picard MarkDuplicates "

# Set output directory
outDir = paste0(projectDir, "data/")


# Helper functions
is_single_end = function(fastq) { fq = get0(fastq, envir = .GlobalEnv); return(!is.null(fq) && !is.na(fq) && fq != "") }
is_paired_end = function(R1, R2) { return(is_single_end(R1) && is_single_end(R2)) }

# Name of STAR's output
starBamFile = "Aligned.out.bam"

for (i in 1:nrow(files)) {
  
  if (!files$include[i]) { next }
  
  time_start = Sys.time()
  print(paste(files$sample[i], time_start))
  
  print("Copy files to data dir") # Record where the data was copied from first
  sampleDir = paste0(outDir, files$sample[i], "/")
  if (!file.exists(paste0(sampleDir, files$fqMate1[i])) || !file.exists(paste0(sampleDir, files$fqMate2[i]))) {
    files$fqMate1_input[i] = paste0(files$dir[i], stringr::str_split_1(files$fqMate1[i], "\\|"), collapse = "|")
    files$fqMate2_input[i] = paste0(files$dir[i], stringr::str_split_1(files$fqMate2[i], "\\|"), collapse = "|")
  }
  mvFiles(files$fqMate1[i], files$dir[i], sampleDir)
  files$dir[i] = mvFiles(files$fqMate2[i], files$dir[i], sampleDir)
  
  print("Format and concatenate fastq files")
  files$fqMate1[i] = formatFastq(files$fqMate1[i], files$dir[i], paste0(files$sample[i], "_1"), threads = threads)
  files$fqMate2[i] = formatFastq(files$fqMate2[i], files$dir[i], paste0(files$sample[i], "_2"), threads = threads)
  
  print("Fastqc")
  system(paste0("fastqc ", files$dir[i],files$fqMate1[i], " -o ", files$dir[i]))
  system(paste0("fastqc ", files$dir[i],files$fqMate2[i], " -o ", files$dir[i]))
  
  print("Cut 3 prime adapter sequences") # These files are only temporary and are removed after mapping
  pelist = fqPECutadapt(files$fqMate1[i], files$fqMate2[i], files$dir[i], seq_platform, threads = threads)
  fqMate1_trimmed = pelist[[1]]
  fqMate2_trimmed = pelist[[2]]
  
  genomeDir = "/home/Apps/SNPsplit/redo/STAR_mm39_N-masked/"
  
  print("STAR mapping")
  star_command = paste("STAR --runThreadN 32 --outSAMtype BAM Unsorted --alignEndsType EndToEnd --outSAMattributes NH HI NM MD --genomeDir", genomeDir, "--outFileNamePrefix",files$dir[i],"--readFilesIn", paste0(files$dir[i], fqMate1_trimmed, " ", files$dir[i], fqMate2_trimmed))
  
  print("SNPsplit-ing")
  bamFile = paste0(files$sample[i],".Aligned.out.bam")
  snpsplit_command <- paste("SNPsplit --snp_file /home/Apps/SNPsplit/refFiles/mm39/all_SNPs_SPRET_EiJ_GRCm39_header.txt", paste0(files$dir[i],bamFile),"--paired")
  
  # For B6 genome
  print("Sort B6 bam")
  b6bamFile = paste0(files$sample[i],".Aligned.out.genome1.bam")
  b6bamSorted = sortBam(paste0(files$dir[i], b6bamFile),
                             bamSortFile = paste0(files$dir[i], files$sample[i], ".genome1.sort"),
                             delBam = F, threads = threads, mem = sortMem)
  
  print("Mark duplicates")
  files$B6bamFile[i] = markDups(paste0(files$dir[i], b6bamSorted), picardCmd, delBam = F)
  
  print("Indexing")
  indexCmd = paste("samtools index", files$B6bamFile[i])
  
  files$B6bamFile[i] = gsub("genome1", "B6", files$B6bamFile[i])
  
  # For Spret genome
  print("Sort Spret bam")
  spretbamFile = paste0(files$sample[i],".Aligned.out.genome2.bam")
  spretbamSorted = sortBam(paste0(files$dir[i], spretbamFile),
                        bamSortFile = paste0(files$dir[i], files$sample[i], ".genome2.sort"),
                        delBam = F, threads = threads, mem = sortMem)
  
  print("Mark duplicates")
  files$spretbamFile[i] = markDups(paste0(files$dir[i], spretbamSorted), picardCmd, delBam = F)
  
  print("Indexing")
  indexCmd = paste("samtools index", files$spretbamFile[i])
  
  files$spretbamFile[i] = gsub("genome2", "spret", files$spretbamFile[i])
  
  # For All reads
  print("Sort allele flagged bam")
  allbamFile = paste0(files$sample[i],".Aligned.out.allele_flagged.bam")
  allbamSorted = sortBam(paste0(files$dir[i], allbamFile),
                           bamSortFile = paste0(files$dir[i], files$sample[i], ".all.sort"),
                           delBam = F, threads = threads, mem = sortMem)
  
  print("Mark duplicates")
  files$allbamFile[i] = markDups(paste0(files$dir[i], allbamSorted), picardCmd, delBam = F)
  
  print("Indexing")
  indexCmd = paste("samtools index", files$allbamFile[i])
  
  files$allbamFile[i] = gsub("all", "all", files$spretbamFile[i])
  
  print("Get mapping stats")
  cts = getBamCts(paste0(files$dir[i], files$allFile[i]))
  
  
  files$unmapped.reads[i] = cts[1]
  files$mapped.reads[i] = cts[2]
  files$unique.reads[i] = cts[3]
  files$paired.reads[i] = cts[4]
  
  # Calculate time to run this sample
  files$timeTaken[i] = difftime(Sys.time(), time_start, units = "hours")
  print(paste("Time taken for sample", files$sample[i], "(hrs):", files$timeTaken[i]))
  
  
  # Update manifest file and note that this sample finished mapping by setting 'include' to FALSE
  files$include[i] = FALSE
  write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F)
}

# Update manifest file
files$include = TRUE
write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F)
####################
# Plot mapping stats

plot_mapping_stats(files)

  
  
  
  
