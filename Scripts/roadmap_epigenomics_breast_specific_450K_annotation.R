#######################################
# Breast (myoepithelial & HMEC) genomic feature overlap w/ 450K array
# Author: Owen Wilkins
#######################################

rm(list = ls())
library(GenomicRanges)
library(rtracklayer)
#BiocInstaller::biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
setwd("/Users/Owen 1/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load annotation data and process to desired format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(Manifest)
data(Locations)
annot <- Manifest
loc <- Locations
# check CGs in same order before combining
all(rownames(loc)==rownames(annot))
# combine 
annot_2 <- cbind(annot, loc)

# Generate new variables for location 
annot_2$hg19_start = annot_2$pos
annot_2$hg19_end = annot_2$pos + 1
# Make the 'CHR' integer variable a character
annot_2$chr = as.character(annot_2$chr)
names(annot_2) 
annot_sub = annot_2[, c("Name", "chr", "hg19_start", "hg19_end")]

# Create a 'GRanges' object from the Illumina annotation file 
Illumina_gr <- makeGRangesFromDataFrame(annot_sub, keep.extra.columns=TRUE, ignore.strand=TRUE, seqnames.field = "chr", start.field="hg19_start", end.field="hg19_end")
# Print the GRange object to see what it looks like
names(Illumina_gr) <- annot_sub$Name
Illumina_gr

# clean workspace
rm(annot, loc, Manifest, Locations, annot_sub)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add genomic feature overlap data of 450K CpGs w/ variuous marks available in Roadmap ChIP-seq data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######## E027 - Breast myoepithelial cells ######## 

#Note on consolidated ChIP-seq data:
#Signal processing engine of MACSv2.0.10 peak caller to generate genome-wide signal coverage tracks
#Negative log10 of the Poisson p-value of ChIP-seq or DNase counts relative to expected background counts
#These signal confidence scores provides a measure of statistical significan of the observed enrichment.

# get ChIP-seq data for available chromatin marks: save paths to data online
#E027_paths <- list("http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K4me1.pval.signal.bigwig", 
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K4me3.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K9ac.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K9me3.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K27me3.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E027-H3K36me3.pval.signal.bigwig")

# import BIGWIG files from consolidated epigenomes from ROADMAP website 
#E027_data <- lapply(E027_paths, import.bw)

# save data so we don't have to download again 
#save(E027_data, file = "/Users/Owen/Dropbox (Christensen Lab)/NDRI_Breast_5hmC_update/05.Enrichment_Analyses/roadmap_breast_myoepith_H3_ChIP_seq.Rdata")
load("03.Enrichment_Analyses/Roadmap_Epigenomics/roadmap_breast_myoepith_H3_ChIP_seq.Rdata")
names(E027_data) <- c("H3K4me1", "H3K4me3", "H3K9ac", "H3K9me3", "H3K27me3", "H3K36me3")

# write a function that will return binary vector of overlap between Granges object and annot rows  
add_annot <- function(mark_object){
  # index for sites w/ -log10P > 2 as roadmap advise this provides good signal/noice separation 
  mark_bs <- mark_object[score(mark_object)>=2, ]
  # find the overlapping regions between the high confidence peaks and the CpGs on the 450k array 
  overlaps <- findOverlaps(mark_bs, Illumina_gr)
  # get the indicies of the overlapping sites in the 450K annotation file 
  indicies <- subjectHits(overlaps)
  # add dummy variable to annotation file for CpGs overlapping with this mark 
  mark_key <- rep(0, nrow(annot_2))
  mark_key[indicies] <- 1
  mark_key
}
# apply fucntion to annotate annot_2
annot_2$Br_myo_H3K4me1 <- add_annot(E027_data$H3K4me1)
annot_2$Br_myo_H3K4me3 <- add_annot(E027_data$H3K4me3)
annot_2$Br_myo_H3K9ac <- add_annot(E027_data$H3K9ac)
annot_2$Br_myo_H3K9me3 <- add_annot(E027_data$H3K9me3)
annot_2$Br_myo_H3K27me3 <- add_annot(E027_data$H3K27me3)
annot_2$Br_myo_H3K36me3 <- add_annot(E027_data$H3K36me3)

# check that it worked by printing count tables 
colnames(annot_2)[14]
apply(annot_2[,14:19], 2, table)

######## E119 - Human Mammary Epithelial cells (HMEC######## 

# E119 = numeric epigenome identifier for consolidated data from HMEC
#E119_paths <- list("http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-DNase.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H2A.Z.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K4me1.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K4me2.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K4me3.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K9ac.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K9me3.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K27ac.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K27me3.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K36me3.pval.signal.bigwig", 
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H3K79me2.pval.signal.bigwig",
#"http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/E119-H4K20me1.pval.signal.bigwig")

#E119_data <- lapply(E119_paths, import.bw)
#save(E119_data, file = "03.Enrichment_Analyses/roadmap_HMEC_ChIP_seq.Rdata")
load("03.Enrichment_Analyses/Roadmap_Epigenomics/roadmap_HMEC_ChIP_seq.Rdata")
names(E119_data) <- c("H2A.Z", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9ac", "H3K9me3", 
                      "H3K27ac", "H3K27me3", "H3K36me3", "H3K79me2", "H3K20me1")

# apply fucntion to annotate annot_2
annot_2$HMEC_H2A.Z <- add_annot(E119_data$H2A.Z)
annot_2$HMEC_H3K4me1 <- add_annot(E119_data$H3K4me1)
annot_2$HMEC_H3K4me2 <- add_annot(E119_data$H3K4me2)
annot_2$HMEC_H3K4me3 <- add_annot(E119_data$H3K4me3)
annot_2$HMEC_H3K9ac <- add_annot(E119_data$H3K9ac)
annot_2$HMEC_H3K9me3 <- add_annot(E119_data$H3K9me3)
annot_2$HMEC_H3K27ac <- add_annot(E119_data$H3K27ac)
annot_2$HMEC_H3K27me3 <- add_annot(E119_data$H3K27me3)
annot_2$HMEC_H3K36me3 <- add_annot(E119_data$H3K36me3)
annot_2$HMEC_H3K79me2 <- add_annot(E119_data$H3K79me2)
annot_2$HMEC_H3K20me1 <- add_annot(E119_data$H3K20me1)

# check that it worked by printing count tables 
colnames(annot_2)[20]
apply(annot_2[,20:30], 2, table)

# save annotated file 
saveRDS(annot_2, "Misc/Illumina-Human-Methylation-450kilmn12-hg19.annotated.rds")
write.csv(annot_2, "Misc/Illumina-Human-Methylation-450kilmn12-hg19.annotated.csv")

