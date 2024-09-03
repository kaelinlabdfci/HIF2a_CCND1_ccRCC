# Loading libraries
library(purrr)
library(dplyr)
library(eulerr)
library(DESeq2)
library(ggplot2)
library(viridis)
library(magrittr)
library(openxlsx)
library(ggsignif)
library("biomaRt")
library(tidyverse)
require(ggVennDiagram)
library(org.Hs.eg.db)




# Merge files function
multmerge = function(mypath){
  filenames=list.files(path=mypath, full.names=TRUE,pattern = "\\.txt$")
  datalist = lapply(filenames, function(x){read.csv(file=x,header=T,skip = 1,sep = "\t")})
  Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}

##################################################################################################################
##################################################################################################################

## Create read count table
OSRC2_TUHR4TKB_24HR_48HR <- multmerge("RNA-seq data/30-601261837/HISAT2_featurecounts/counts")
OSRC2_TUHR4TKB_24HR_48HR_subsetted <- OSRC2_TUHR4TKB_24HR_48HR %>% dplyr::select(Geneid,
                                                                          contains("X.omics.odcf.analysis."))
colnames(OSRC2_TUHR4TKB_24HR_48HR_subsetted) <- gsub("X.omics.odcf.analysis.OE0497_projects.pedHGG_ChIP.analysis_d132r.30.601261837.|.sorted.bam","",colnames(OSRC2_TUHR4TKB_24HR_48HR_subsetted))
colnames(OSRC2_TUHR4TKB_24HR_48HR_subsetted) <- c("Geneid","OSRC2-DMSO-1","4TKB-DMSO-1",
                                                  "4TKB-DMSO-2","4TKB-DMSO-3",
                                                  "4TKB-PT24h-1","4TKB-PT24h-2",
                                                  "4TKB-PT24h-3","4TKB-PT48h-1",
                                                  "4TKB-PT48h-2","4TKB-PT48h-3",
                                                  "OSRC2-DMSO-2","OSRC2-DMSO-3",
                                                  "OSRC2-PT24h-1","OSRC2-PT24h-2",
                                                  "OSRC2-PT24h-3","OSRC2-PT48h-1",
                                                  "OSRC2-PT48h-2","OSRC2-PT48h-3"
                                                  )
OSRC2_TUHR4TKB_24HR_48HR_subsetted_summarise <- OSRC2_TUHR4TKB_24HR_48HR_subsetted %>%
  mutate(Geneid= sub("[.][0-9]*","",Geneid))%>%
  group_by(Geneid) %>%
  summarise_all(mean, na.rm = TRUE) %>% as.data.frame()

# gene IDs should be stored as row.names
row.names(OSRC2_TUHR4TKB_24HR_48HR_subsetted_summarise) <- OSRC2_TUHR4TKB_24HR_48HR_subsetted_summarise$Geneid
OSRC2_TUHR4TKB_24HR_48HR_subsetted_summarise$Geneid <- NULL

# Create a sample information dataframe for DESeq2 analysis
sample_info <- DataFrame(condition = gsub("-[0-9]+", "", names(OSRC2_TUHR4TKB_24HR_48HR_subsetted_summarise)), row.names = names(OSRC2_TUHR4TKB_24HR_48HR_subsetted_summarise) )

# DESeq2 object construction
DESeq.ds <- DESeqDataSetFromMatrix(countData = OSRC2_TUHR4TKB_24HR_48HR_subsetted_summarise, 
                                   colData = sample_info,design = ~ condition)

colSums(counts(DESeq.ds)) %>% barplot

# Only keep genes that are expressed atleast in 1 sample
keep_genes <- rowSums(counts(DESeq.ds)) > 0 
DESeq.ds <- DESeq.ds[ keep_genes, ] 

#Now that we have the data, we can start using DESeq’s functions, e.g. estimateSizeFactors() for sequencing depth normalization.
plot(sizeFactors(DESeq.ds), colSums(counts(DESeq.ds)))

DESeq.ds <- estimateSizeFactors(DESeq.ds) 
sizeFactors(DESeq.ds)


# Checking whether the normalization helped adjust global differences between the samples.
# SF normalized read counts only

# setting up the plotting layout
par(mfrow=c(1,2))
counts.sf_normalized <- log2(counts(DESeq.ds, normalized=TRUE)+1)
# adding the boxplots
boxplot(counts.sf_normalized, main = "SF normalized") 
boxplot(counts(DESeq.ds), main = "read counts only")

write.csv(counts(DESeq.ds, normalize= TRUE), "Processed_data/Normalised_counts_OSRC2_TUHR_24_48_PT_DMSO_DESeq2.csv")
write.csv(counts(DESeq.ds), "Processed_data/Non-normalised_counts_OSRC2_TUHR_24_48_PT_DMSO.csv")

par(mfrow=c(1,2)) # to plot the two box plots next to each other 
boxplot(log2(counts(DESeq.ds)), notch=TRUE,
main = "Non-normalized read counts\n(log-transformed)",
ylab="read counts")
boxplot(log2(counts(DESeq.ds, normalize= TRUE)), notch=TRUE,
        main = "Size-factor-normalized read counts\n(log-transformed)", ylab="read counts")

dds <- DESeq(DESeq.ds)

# PT vs DMSO comparisons
### RNAseq PT vs DMSO OSRC2 24hr 
res.lfc.PT.DMSO.OSRC2.24 <- results(dds,contrast = c("condition","OSRC2-PT24h", "OSRC2-DMSO")) %>% as.data.frame() %>% 
  na.omit()
ens.str <- substr(rownames(res.lfc.PT.DMSO.OSRC2.24), 1, 15)
res.lfc.PT.DMSO.OSRC2.24$symbol <- mapIds(org.Hs.eg.db,
                                  keys=ens.str,
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")
write.csv(res.lfc.PT.DMSO.OSRC2.24%>%as.data.frame(),"Processed_data/RNAseq_OSRC2_24HR_PT_DMSO.csv")

### RNAseq PT vs DMSO OSRC2 48hr 
res.lfc.PT.DMSO.OSRC2.48 <- results(dds,contrast = c("condition","OSRC2-PT48h", "OSRC2-DMSO")) %>% as.data.frame() %>% 
  na.omit()
ens.str <- substr(rownames(res.lfc.PT.DMSO.OSRC2.48), 1, 15)
res.lfc.PT.DMSO.OSRC2.48$symbol <- mapIds(org.Hs.eg.db,
                                          keys=ens.str,
                                          column="SYMBOL",
                                          keytype="ENSEMBL",
                                          multiVals="first")
write.csv(res.lfc.PT.DMSO.OSRC2.48%>%as.data.frame(),"Processed_data/RNAseq_OSRC2_48HR_PT_DMSO.csv")

### RNAseq PT vs DMSO TUHR4TKB 24hr 
res.lfc.PT.DMSO.TUHR.24 <- results(dds,contrast = c("condition","4TKB-PT24h", "4TKB-DMSO")) %>% as.data.frame() %>% 
  na.omit()
ens.str <- substr(rownames(res.lfc.PT.DMSO.TUHR.24), 1, 15)
res.lfc.PT.DMSO.TUHR.24$symbol <- mapIds(org.Hs.eg.db,
                                          keys=ens.str,
                                          column="SYMBOL",
                                          keytype="ENSEMBL",
                                          multiVals="first")
write.csv(res.lfc.PT.DMSO.TUHR.24%>%as.data.frame(),"Processed_data/RNAseq_TUHR4TKB_24HR_PT_DMSO.csv")

### RNAseq PT vs DMSO TUHR4TKB 48hr 
res.lfc.PT.DMSO.TUHR.48 <- results(dds,contrast = c("condition","4TKB-PT48h", "4TKB-DMSO")) %>% as.data.frame() %>% 
  na.omit()
ens.str <- substr(rownames(res.lfc.PT.DMSO.TUHR.48), 1, 15)
res.lfc.PT.DMSO.TUHR.48$symbol <- mapIds(org.Hs.eg.db,
                                         keys=ens.str,
                                         column="SYMBOL",
                                         keytype="ENSEMBL",
                                         multiVals="first")
write.csv(res.lfc.PT.DMSO.TUHR.48%>%as.data.frame(),"Processed_data/RNAseq_TUHR4TKB_48HR_PT_DMSO.csv")
