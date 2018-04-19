######################
# Add breast-specific super enhancer and 
# regulatory regions to annotation file
#
# Author: Alexander Titus
# Created: 04/15/2018
######################

######################
     # Set up working environment
     #install.packages('CREAM')
     library(CREAM)
     library(data.table)
     working.dir = 'C:/Users/atitus/github/shared_450K_annotation'
     setwd(working.dir)

     
######################
     # Read in annotation file
     file.name = 'Illumina-Human-Methylation-450kilmn12-hg19.annotated.rds'
     anno450k = data.frame(readRDS(file.name))
     
     
######################
     # Download super enhancer data from dbSuper
     hmec.url = 'http://asntech.org/dbsuper/data/bed/hg19/HMEC.bed'
     mcf7.url = 'http://asntech.org/dbsuper/data/bed/hg19/MCF-7.bed'
     hmec.file = 'data/HMEC_dbSuper.bed'
     mcf7.file = 'data/MCF7_dbSuper.bed'
     download.file(hmec.url, destfile = hmec.file)
     download.file(mcf7.url, destfile = mcf7.file)
     
     
######################
     # Add super enhancer data to the annotation file
     hmec.anno = data.frame(read.table(hmec.file, header = FALSE, sep="\t", stringsAsFactors=FALSE, quote=""))
     colnames(hmec.anno) = c('chr', 'start', 'end', 'ID', 'count')
     
     anno450k$HMEC_SuperEnhancer = rep(0, nrow(anno450k))
     anno450k$HMEC_SuperEnhancer = ifelse(anno450k$pos >= hmec.anno$start & anno450k$pos <= hmec.anno$end, 1, 0)

     mcf7.anno = data.frame(read.table(mcf7.file, header = FALSE, sep="\t", stringsAsFactors=FALSE, quote=""))
     colnames(mcf7.anno) = c('chr', 'start', 'end', 'ID', 'count')
     
     anno450k$MCF7_SuperEnhancer = rep(0, nrow(anno450k))
     anno450k$MCF7_SuperEnhancer = ifelse(anno450k$pos >= mcf7.anno$start & anno450k$pos <= mcf7.anno$end, 1, 0)
     
     saveRDS(anno450k, "Illumina-Human-Methylation-450kilmn12-hg19.annotated.rds")
     csv.file = 'E:/Dropbox (Christensen Lab)/Shared 450K annotation data/Illumina-Human-Methylation-450kilmn12-hg19.annotated.csv'
     fwrite(anno450k, csv.file)
     
     
######################
     # Calculate CORES
     # https://github.com/bhklab/CREAM
     hmec.chromeaccess = 'C:/Users/atitus/Downloads/UCSF-UBC.Breast_Luminal_Epithelial_Cells.Input.RM080.bed/UCSF-UBC.Breast_Luminal_Epithelial_Cells.Input.RM080.bed'
     hmec.COREfile = 'data/HMEC_CORE.csv'
     mcf7.COREfile = 'data/MCF7_CORE.csv'
     
     CREAM( in_path = hmec.chromeaccess, 
            out_path = hmec.COREfile, MinLength = 1000, peakNumMin = 2)
     