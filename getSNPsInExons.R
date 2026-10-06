library("GenomicRanges")
library("GenomicFeatures")
library(VariantAnnotation)

# Read in bistools; select_ref_txdb, select_ref_org, filter_detected
source("/home/Apps/bitbucket/bistools/ESB_bisTools.R") # Magic

genome = "mm39"
tx = select_ref_txdb(genome)
org = select_ref_org(genome)

# Set up exon and gene tables for annotation of reads
exons = exonsBy(tx, by = "gene")
ex = as.data.frame(exons)
colnames(ex) = c("group", "ENTREZID", "chr", "start", "end", "width", "strand", "exon_id", "exon_name")
ex$chr = gsub("chr","",ex$chr)

# get gene names
geneSym = AnnotationDbi::select(org, columns = "SYMBOL", keytype = "ENTREZID", keys = names(exons))

# merge the exon info with gene symobols
ex_merged <- merge(ex, geneSym, by.x = "ENTREZID", by.y = "ENTREZID", all.x = TRUE)
ex_merged$chr <- sub("^chr", "", ex_merged$chr)

gr <- makeGRangesFromDataFrame(ex_merged)

# read in variant file
# vcf_file <- "/home/boss_lab/Apps/genomes/species/spret-b6/SPRET_EiJ.variants/SPRET_EiJ.mgp.v5.snps.dbSNP142.vcf"
# vr <- readVcfAsVRanges(vcf_file)
variants_file = "/home/Apps/SNPsplit/refFiles/mm39/all_SNPs_SPRET_EiJ_GRCm39_header.txt"
variants <- read.table(variants_file, sep = "\t", header = TRUE, as.is = TRUE)
variants$start <- variants$Position
variants$end <- variants$Position
variants <- variants[, c("Chromosome", "start", "end", "Strand", "Ref.SNP", "SNP.ID")]
variants$Strand = "*"

vr2 = makeGRangesFromDataFrame(variants, keep.extra.columns = TRUE)

# Which variants overlap with the exonic regions
# find overlaps
over2 = findOverlaps(gr, vr2)

# Extract the overlapping regions
snp_exon_pairs <- cbind(as.data.frame(gr[queryHits(over2)]),
                        as.data.frame(vr2[subjectHits(over2)]))

# select columns of interest
fin = snp_exon_pairs[,c(1,2,3,6,7,11,12)]
colnames(fin) = c("exon_chr", "exon_start", "exon_end", "snp_chr", "snp_pos", "ref/alt", "snpID")

# fix columns and add unique identifiers
fin$exon_chr = paste0("chr", fin$exon_chr)
fin$snp_chr = paste0("chr", fin$snp_chr)
fin$exon_id = paste0(fin$exon_chr, "_", fin$exon_start,"_", fin$exon_end)
# fin$snp_id = paste0(fin$snp_chr, "_", fin$snp_pos,"_", fin$snp_pos)

# add id column to exon for merging with final df
ex_merged$chr = paste0("chr", ex_merged$chr)
ex_merged$exon_id2 = paste0(ex_merged$chr, "_", ex_merged$start,"_", ex_merged$end)

# merge to get entrezid and gene symbols
fin$ENTREZID <- ex_merged$ENTREZID[match(fin$exon_id, ex_merged$exon_id2)]
fin$SYMBOL <- ex_merged$SYMBOL[match(fin$exon_id, ex_merged$exon_id2)]

setwd("/home/mgupta/snpsplit/MZB_FoB_2023.11/RNA/counts/")
write.table(fin, "SNPs_present_in_exons.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# prep for liftover
# exon table
exon = fin[c(1,2,3,8)]
write.table(exon, "exons_containing_snps.bed", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# snp table
snp = fin[c(4,5,5,9)]
write.table(snp, "snps_in_genes.bed", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)





