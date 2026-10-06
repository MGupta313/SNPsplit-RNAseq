# Binomial Testing for MZB and FoB
setwd("/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/BinomialTest")

# Separate MZB and FoB samples
manifest = read.table("/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/RNAseq.sample.manifest.txt", sep ="\t", header = TRUE)
MZB_samples = manifest$sample[manifest$group == "MZ"]
Fob_samples = manifest$sample[manifest$group == "FO"]

########### Combining all MZB samples ########################
# Read in read counts at gene
geneDir = "/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/readCountsAtGene/"

# List of file names
file_names <- paste0(MZB_samples, "_readCountsAtGene.txt")

# Create an empty list to store the data frames
dfs <- list()

# Loop through each file
for (file_name in file_names) {
  
  # Read the data from the file
  df <- read.table(paste0(geneDir, file_name), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract the sample name (SCR001, SCR002, etc.) from the file name
  sample_name <- gsub("_readCountsAtGene.txt", "", file_name)
  
  # Add the sample name as a prefix to the relevant columns
  colnames(df)[colnames(df) == "B6_reads"] <- paste0(sample_name, "_B6_reads")
  colnames(df)[colnames(df) == "Spret_reads"] <- paste0(sample_name, "_Spret_Reads")
  
  # Append the data frame to the list
  dfs[[file_name]] <- df
}

# Merge all the data frames by ENTREZID and SYMBOL
merged_df <- dfs[[1]]  # Start with the first data frame
for (i in 2:length(dfs)) {
  merged_df <- merge(merged_df, dfs[[i]], by = c("ENTREZID", "SYMBOL"))
}

# Reorder columns to have all B6_reads columns first, followed by Spret_reads columns
b6_columns <- grep("B6_reads", colnames(merged_df), value = TRUE)
spret_columns <- grep("Spret_Reads", colnames(merged_df), value = TRUE)

# Reorder the columns
ordered_columns <- c("ENTREZID", "SYMBOL", b6_columns, spret_columns)

# Reorder the dataframe
merged_df <- merged_df[, ordered_columns]

# View the reordered dataframe
head(merged_df)

MZB_merged_df = merged_df

########### Combining all Fob samples ########################
# Read in read counts at gene
geneDir = "/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/readCountsAtGene/"

# List of file names
file_names <- paste0(Fob_samples, "_readCountsAtGene.txt")

# Create an empty list to store the data frames
dfs <- list()

# Loop through each file
for (file_name in file_names) {
  
  # Read the data from the file
  df <- read.table(paste0(geneDir, file_name), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract the sample name (SCR001, SCR002, etc.) from the file name
  sample_name <- gsub("_readCountsAtGene.txt", "", file_name)
  
  # Add the sample name as a prefix to the relevant columns
  colnames(df)[colnames(df) == "B6_reads"] <- paste0(sample_name, "_B6_reads")
  colnames(df)[colnames(df) == "Spret_reads"] <- paste0(sample_name, "_Spret_Reads")
  
  # Append the data frame to the list
  dfs[[file_name]] <- df
}

# Merge all the data frames by ENTREZID and SYMBOL
merged_df <- dfs[[1]]  # Start with the first data frame
for (i in 2:length(dfs)) {
  merged_df <- merge(merged_df, dfs[[i]], by = c("ENTREZID", "SYMBOL"))
}

# Reorder columns to have all B6_reads columns first, followed by Spret_reads columns
b6_columns <- grep("B6_reads", colnames(merged_df), value = TRUE)
spret_columns <- grep("Spret_Reads", colnames(merged_df), value = TRUE)

# Reorder the columns
ordered_columns <- c("ENTREZID", "SYMBOL", b6_columns, spret_columns)

# Reorder the dataframe
merged_df <- merged_df[, ordered_columns]

# View the reordered dataframe
head(merged_df)

Fob_merged_df = merged_df

########## binomial on MZB ##########
genes = MZB_merged_df[c(1,2)]
mzb = MZB_merged_df[-c(1,2)]

# Create group labels for each sample
group <- c(rep("B6", 4), rep("Spret", 4))

# Perform Two-group beta-binomial test
out <- countdata::bb.test(mzb, 
                          colSums(mzb), 
                          group)

d.norm <- countdata::normalize(mzb)
colnames(d.norm) = paste0("Norm.", colnames(d.norm))


x = cbind(genes,mzb, d.norm, 
          fc = countdata::fold.change(d.norm[, 1:4], d.norm[, 5:8]),
          pval = out$p.value,
          pval.BH = p.adjust(out$p.value, method = "BH"))



write.table(x, file = "MZB_B6.v.Spret_BinomialTest.txt", row.names = FALSE, sep = "\t")


########## binomial on Fob ##########
genes = Fob_merged_df[c(1,2)]
fob = Fob_merged_df[-c(1,2)]

# Create group labels for each sample
group <- c(rep("B6", 4), rep("Spret", 4))

# Perform Two-group beta-binomial test
out <- countdata::bb.test(fob, 
                          colSums(fob), 
                          group)

d.norm <- countdata::normalize(fob)
colnames(d.norm) = paste0("Norm.", colnames(d.norm))


x = cbind(genes,fob, d.norm, 
          fc = countdata::fold.change(d.norm[, 1:4], d.norm[, 5:8]),
          pval = out$p.value,
          pval.BH = p.adjust(out$p.value, method = "BH"))



write.table(x, file = "FoB_B6.v.Spret_BinomialTest.txt", row.names = FALSE, sep = "\t", col.names = TRUE)
