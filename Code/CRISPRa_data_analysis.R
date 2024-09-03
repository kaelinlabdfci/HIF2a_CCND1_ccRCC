## Load libraries
library(purrr)
library(dplyr)
library(eulerr)
library(ggplot2)
library(viridis)
library(magrittr)
library(openxlsx)
library(ggsignif)
library(tidyverse)
require(ggVennDiagram)

## Limma comparison for CRISPRa
run_lm_stats_limma <- function(mat, vec, covars = NULL, weights = NULL, target_type = 'Gene') {
  require(limma)
  require(magrittr)
  require(plyr)
  
  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2) #only two group comparisons implemented so far
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],,drop=F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],,drop=F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  } else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata,]))
  }
  #there must be more than one unique value of the independent variable
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  #if using covariates add them as additional predictors to the model
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata,, drop = FALSE]
    combined[['pred']] <- pred
    form <- as.formula(paste('~', paste0(colnames(combined), collapse = ' + ')))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  } else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata,])
    } else{
      weights <- weights[udata]
    }
  }
  fit <- limma::lmFit(t(mat[udata,]), design, weights = weights, method = 'robust')
  fit <- limma::eBayes(fit)
  targ_coef <- grep('pred', colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf)
  ## logFC column gives the value of the contrast. Usually this represents a log2-fold change
  ## between two or more experimental conditions
  
  ## AveExpr column gives the average log2-expression level for that gene across all the arrays
  ## and channels in the experiment
  
  ## t is the moderated t-statistic
  
  ## Column P.Value is the associated p-value and adj.P.Value is the p-value adjusted for multiple 
  ## testing. Benjamini and Hochberg’s method to control the false discovery rate. 
  ## The adjusted values are often called q-values if the intention is to control or estimate the
  ## false discovery rate.If all genes with q-value below a threshold, say 0.05, are selected as 
  ## differentially expressed, then the expected proportion of false discoveries in the selected 
  ## group is controlled to be less than the threshold value, in this case 5%.
  
  ## B-statistic (lods or B) is the log-odds that the gene is differentially expressed. Suppose 
  ## for example that B = 1.5. The odds of differential expression is exp(1.5)=4.48, i.e, about
  ## four and a half to one. The probability that the gene is differentially expressed is 4.48/(1+4.48)=0.82,
  ## i.e., the probability is about 82% that this gene is differentially expressed. A B-statistic of zero corresponds to a 50-50 chance that the gene is differentially expressed. The B-statistic is automatically
  ## adjusted for multiple testing by assuming that 1% of the genes, or some other percentage specified
  ## are expected to be differentially expressed
  
  if (colnames(results)[1] == 'ID') {
    colnames(results)[1] <- target_type
  } else {
    results %<>% rownames_to_column(var = target_type)
  }
  results$min_samples <- min_samples[results[[target_type]]]
  
  results %<>% set_colnames(revalue(colnames(.), c('logFC' = 'EffectSize', 'AveExpr' = 'Avg', 't' = 't_stat', 'B' = 'log_odds',
                                                   'P.Value' = 'p.value', 'adj.P.Val' = 'q.value', 'min_samples' = 'min_samples'))) %>% na.omit()
  return(results)
}

##################################################################################################################
##################################################################################################################

## data import - TUHR4TKB
barcode_gene_mapping_setA <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/CP0052_origtarget_20191021.chip", 
                                        sep = "\t", header = T)
barcode_gene_mapping_setB <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/CP0053_origtarget_20191021.chip", 
                                        sep = "\t", header = T)

plate7 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-7/DR_GPP4421_2903468_Shirole_Plate7/lognorm-DR_GPP4421_2903468_Shirole_Plate7.txt",
                     sep = "\t", header = T) %>% inner_join(barcode_gene_mapping_setA, 
                                                            by = c("Construct.Barcode"="Barcode.Sequence")) %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.Day.0_TUHR4TKB, Rep.2.Day.0_TUHR4TKB, Empty) 

plate10 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-10/DR_GPP4417_2903468_Shirole_Plate10/lognorm-DR_GPP4417_2903468_Shirole_Plate10.txt",
                      sep = "\t", header = T) %>% inner_join(barcode_gene_mapping_setB, 
                                                             by = c("Construct.Barcode"="Barcode.Sequence")) %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.Day.0_TUHR4TKB, Rep.2.Day.0_TUHR4TKB, Empty) 

plate_8 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-8/DR_GPP4422_2903468_Shirole_Plate08/lognorm-DR_GPP4422_2903468_Shirole_Plate08.txt",
                      sep = "\t", header = T) %>% select(-Empty) %>% full_join(read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-9/DR_GPP4423_2903468_Shirole_Plate09/lognorm-DR_GPP4423_2903468_Shirole_Plate09.txt",
                                                                                          sep = "\t", header = T)) %>% 
  inner_join(barcode_gene_mapping_setA, by = c("Construct.Barcode"="Barcode.Sequence")) 
platedmso <- plate_8 %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.DMSO.Day.20_TUHR4TKB, Rep.2.DMSO.Day.20_TUHR4TKB) 
plate_s <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-11/JD_GPP4418_2903468_Shirole_Plate11/lognorm-JD_GPP4418_2903468_Shirole_Plate11.txt",
                      sep = "\t", header = T) %>% select(-Empty) %>% full_join(read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-12/JD_GPP4419_2903468_Shirole_Plate12/lognorm-JD_GPP4419_2903468_Shirole_Plate12.txt",
                                                                                          sep = "\t", header = T)) %>% inner_join(barcode_gene_mapping_setB, 
                                                                                                                                  by = c("Construct.Barcode"="Barcode.Sequence")) 
plate11 <- plate_s %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.DMSO.Day.20_TUHR4TKB, Rep.2.DMSO.Day.20_TUHR4TKB) 


platept <- plate_8 %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.PT.Day.20_TUHR4TKB, Rep.2.PT.Day.20_TUHR4TKB) 

plate12 <- plate_s %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.PT.Day.20_TUHR4TKB, Rep.2.PT.Day.20_TUHR4TKB) 


sgRNA_abundance_DMSOd20_comb <- rbind(platedmso %>% select(geneSymbol, Construct.IDs, Rep.1.DMSO.Day.20_TUHR4TKB,Rep.2.DMSO.Day.20_TUHR4TKB),
                                      plate11 %>% select(geneSymbol, Construct.IDs, Rep.1.DMSO.Day.20_TUHR4TKB,Rep.2.DMSO.Day.20_TUHR4TKB))

sgRNA_abundance_PTd20_comb <- rbind(platept %>% select(geneSymbol, Construct.IDs, Rep.1.PT.Day.20_TUHR4TKB, Rep.2.PT.Day.20_TUHR4TKB),
                                    plate12 %>% select(geneSymbol, Construct.IDs, Rep.1.PT.Day.20_TUHR4TKB, Rep.2.PT.Day.20_TUHR4TKB))
sgRNA_abundance_d0_comb <- rbind(plate7,plate10) %>% select(-Empty)

genes_not_to_include <- c("MULTIPLE", "POTENTIALLY", "LOC", "LINC", "INACTIVE","NO_SITE","ONE_INTERGENIC_SITE")


### sgRNA abundance DMSO-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_DMSOd20_comb,sgRNA_abundance_d0_comb) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.DMSO.day0.comb = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_comb,sgRNA_abundance_d0_comb) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.day0.comb = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT vs DMSO d20 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_comb,sgRNA_abundance_DMSOd20_comb) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.DMSO.comb = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

write.csv(res.lfc.DMSO.day0.comb,"Processed_data/TUHR4TKB_DMSO_vs_Day0_CRISPRa.csv")
write.csv(res.lfc.PT.day0.comb,"Processed_data/TUHR4TKB_PT_vs_Day0_CRISPRa.csv")
write.csv(res.lfc.PT.DMSO.comb,"Processed_data/TUHR4TKB_PT_vs_DMSO_CRISPRa.csv")


## data import - OSRC2
# setA
plate1 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-1/lognorm-JD_GPP3205_Shirole_Plate1.txt",
                     sep = "\t", header = T) %>% inner_join(barcode_gene_mapping_setA, 
                                                            by = c("Construct.Barcode"="Barcode.Sequence")) %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.Day.0, Rep.2.Day.0, Empty)
plate2 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-2/lognorm-JD_GPP3206_Shirole_Plate2.txt",
                     sep = "\t", header = T) %>% inner_join(barcode_gene_mapping_setA, 
                                                            by = c("Construct.Barcode"="Barcode.Sequence")) %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.DMSO.Day.20, Rep.1.PT.Day.20, Empty) 
plate3 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-3/lognorm-JD_GPP3207_Shirole_Plate3.txt",
                     sep = "\t", header = T) %>% inner_join(barcode_gene_mapping_setA, 
                                                            by = c("Construct.Barcode"="Barcode.Sequence")) %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.2.DMSO.Day.20, Rep.2.PT.Day.20, Empty) 


# setB
plate4 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-4/lognorm-JD_GPP3454_Shirole_20220510_P4.txt",
                     sep = "\t", header = T) %>% inner_join(barcode_gene_mapping_setB, 
                                                            by = c("Construct.Barcode"="Barcode.Sequence")) %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.Day.0, Rep.2.Day.0, Empty) 
plate5 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-5/lognorm-JD_GPP3455_Shirole_20220510_P5.txt",
                     sep = "\t", header = T) %>% inner_join(barcode_gene_mapping_setB, 
                                                            by = c("Construct.Barcode"="Barcode.Sequence")) %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.1.DMSO.Day.20, Rep.1.PT.Day.20, Empty) 
plate6 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-6/lognorm-JD_GPP3456_Shirole_20220510_P6.txt",
                     sep = "\t", header = T) %>% inner_join(barcode_gene_mapping_setB, 
                                                            by = c("Construct.Barcode"="Barcode.Sequence")) %>%
  select(geneSymbol = Annotated.Gene.Symbol, Construct.IDs, Rep.2.DMSO.Day.20, Rep.2.PT.Day.20, Empty) 


sgRNA_abundance_DMSOd20_comb <- full_join(rbind(plate2 %>% select(geneSymbol, Construct.IDs, Rep.1.DMSO.Day.20),
                                                plate5 %>% select(geneSymbol, Construct.IDs, Rep.1.DMSO.Day.20)), 
                                          rbind(plate3 %>% select(geneSymbol, Construct.IDs, Rep.2.DMSO.Day.20),
                                                plate6 %>% select(geneSymbol, Construct.IDs, Rep.2.DMSO.Day.20))
)
sgRNA_abundance_PTd20_comb <- full_join(rbind(plate2 %>% select(geneSymbol, Construct.IDs, Rep.1.PT.Day.20), 
                                              plate5 %>% select(geneSymbol, Construct.IDs, Rep.1.PT.Day.20)),
                                        rbind(plate3 %>% select(geneSymbol, Construct.IDs, Rep.2.PT.Day.20),
                                              plate6 %>% select(geneSymbol, Construct.IDs, Rep.2.PT.Day.20)))
sgRNA_abundance_d0_comb <- rbind(plate1,plate4) %>% select(-Empty)


### sgRNA abundance DMSO-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_DMSOd20_comb,sgRNA_abundance_d0_comb) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.DMSO.day0.comb.OSRC2 = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_comb,sgRNA_abundance_d0_comb) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.day0.comb.OSRC2 = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT vs DMSO d20 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_comb,sgRNA_abundance_DMSOd20_comb) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.DMSO.comb.OSRC2 = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

write.csv(res.lfc.DMSO.day0.comb.OSRC2,"Processed_data/OSRC2_DMSO_vs_Day0_CRISPRa.csv")
write.csv(res.lfc.PT.day0.comb.OSRC2,"Processed_data/OSRC2_PT_vs_Day0_CRISPRa.csv")
write.csv(res.lfc.PT.DMSO.comb.OSRC2,"Processed_data/OSRC2_PT_vs_DMSO_CRISPRa.csv")

## data import - OSRC2 mini library
# /Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-15/lognorm-JD_GPP4898_2903697_Shirole_plate_14.txt
barcode_mapping_miniset <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-15/CP1904_GRCh38_NCBI_CRISPRa_strict_gene_20230726 (2).chip", 
                                      sep = "\t", header = T)
plate15 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-15/lognorm-JD_GPP4898_2903697_Shirole_plate_14.txt",
                     sep = "\t", header = T) %>% inner_join(barcode_mapping_miniset, 
                                                            by = c("Construct.Barcode"="Barcode.Sequence"))

sgRNA_abundance_DMSOd20_miniset <- plate15 %>% select(geneSymbol = Gene.Symbol, Construct.IDs, 
                                                      OSRC2.DMSO.Day.20.Rep.1, OSRC2.DMSO.Day.20.Rep.2,OSRC2.DMSO.Day.20.Rep.3)

sgRNA_abundance_PTd20_miniset <- plate15 %>% select(geneSymbol = Gene.Symbol, Construct.IDs, 
                                                    OSRC2.PT.Day.20.Rep.1, OSRC2.PT.Day.20.Rep.2,OSRC2.PT.Day.20.Rep.3)
sgRNA_abundance_d0_miniset <- plate15 %>% select(geneSymbol = Gene.Symbol, Construct.IDs, 
                                                 OSRC2.Day.0.Rep.1, OSRC2.Day.0.Rep.2,OSRC2.Day.0.Rep.3)


### sgRNA abundance DMSO-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_DMSOd20_miniset,sgRNA_abundance_d0_miniset) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,1,0,0,0))
res.lfc.DMSO.day0.miniset.OSRC2 = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_miniset,sgRNA_abundance_d0_miniset) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,1,0,0,0))
res.lfc.PT.day0.miniset.OSRC2 = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT vs DMSO d20 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_miniset,sgRNA_abundance_DMSOd20_miniset) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,1,0,0,0))
res.lfc.PT.DMSO.miniset.OSRC2 = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

write.csv(res.lfc.DMSO.day0.miniset.OSRC2,"Processed_data/OSRC2_DMSO_vs_Day0_CRISPRa_HIF2a_library.csv")
write.csv(res.lfc.PT.day0.miniset.OSRC2,"Processed_data/OSRC2_PT_vs_Day0_CRISPRa_HIF2a_library.csv")
write.csv(res.lfc.PT.DMSO.miniset.OSRC2,"Processed_data/OSRC2_PT_vs_DMSO_CRISPRa_HIF2a_library.csv")

## data import - TUHR4TKB mini library
barcode_mapping_miniset <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-13/CP1904_GRCh38_NCBI_CRISPRa_strict_gene_20230726.chip", 
                                      sep = "\t", header = T)

plate13 <- read.table("/Users/devishi/Dropbox/Nitin_Chipseq/CRISPRa/Screen Results/Plate-13/JD_GPP4420_2903468_Shirole_Plate13/lognorm-JD_GPP4420_2903468_Shirole_Plate13.txt",
                      sep = "\t", header = T) %>% inner_join(barcode_mapping_miniset, 
                                                             by = c("Construct.Barcode"="Barcode.Sequence")) %>% select(-Gene.ID)

sgRNA_abundance_DMSOd20_miniset <- plate13 %>% select(geneSymbol = Gene.Symbol, Construct.IDs, 
                                                      TUHR.DMSO.Day.20.Rep.1 = X4TKB.DMSO.Day.20.Rep.1, 
                                                      TUHR.DMSO.Day.20.Rep.2 = X4TKB.DMSO.Day.20.Rep.2)

sgRNA_abundance_PTd20_miniset <- plate13 %>% select(geneSymbol = Gene.Symbol, Construct.IDs, 
                                                    TUHR.PT.Day.20.Rep.1 = X4TKB.PT.Day.20.Rep.1, 
                                                      TUHR.PT.Day.20.Rep.2 = X4TKB.PT.Day.20.Rep.2)
sgRNA_abundance_d0_miniset <- plate13 %>% select(geneSymbol = Gene.Symbol, Construct.IDs, 
                                                 TUHR.Day.0.Rep.1 = X4TKB.Day.0.Rep.1, 
                                                 TUHR.Day.0.Rep.2 = X4TKB.Day.0.Rep.2)


### sgRNA abundance DMSO-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_DMSOd20_miniset,sgRNA_abundance_d0_miniset) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.DMSO.day0.miniset.TUHR = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_miniset,sgRNA_abundance_d0_miniset) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.day0.miniset.TUHR = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT vs DMSO d20 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_miniset,sgRNA_abundance_DMSOd20_miniset) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.DMSO.miniset.TUHR = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

write.csv(res.lfc.DMSO.day0.miniset.TUHR,"Processed_data/TUHR4TKB_DMSO_vs_Day0_CRISPRa_HIF2a_library.csv")
write.csv(res.lfc.PT.day0.miniset.TUHR,"Processed_data/TUHR4TKB_PT_vs_Day0_CRISPRa_HIF2a_library.csv")
write.csv(res.lfc.PT.DMSO.miniset.TUHR,"Processed_data/TUHR4TKB_PT_vs_DMSO_CRISPRa_HIF2a_library.csv")

## data import - 786O mini library
sgRNA_abundance_DMSOd20_miniset <- plate13 %>% select(geneSymbol = Gene.Symbol, Construct.IDs, 
                                                      O786.DMSO.Day.20.Rep.1 = X786O.DMSO.Day.20.Rep.1, 
                                                      O786.DMSO.Day.20.Rep.2 = X786O.DMSO.Day.20.Rep.2)

sgRNA_abundance_PTd20_miniset <- plate13 %>% select(geneSymbol = Gene.Symbol, Construct.IDs, 
                                                    O786.PT.Day.20.Rep.1 = X786O.PT.Day.20.Rep.1, 
                                                    O786.PT.Day.20.Rep.2 = X786O.PT.Day.20.Rep.2)
sgRNA_abundance_d0_miniset <- plate13 %>% select(geneSymbol = Gene.Symbol, Construct.IDs, 
                                                 O786.Day.0.Rep.1 = X786O.Day.0.Rep.1, 
                                                 O786.Day.0.Rep.2 = X786O.Day.0.Rep.2)


### sgRNA abundance DMSO-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_DMSOd20_miniset,sgRNA_abundance_d0_miniset) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.DMSO.day0.miniset.786O = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT-d20 vs d0 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_miniset,sgRNA_abundance_d0_miniset) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.day0.miniset.786O = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

### sgRNA abundance PT vs DMSO d20 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_miniset,sgRNA_abundance_DMSOd20_miniset) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.DMSO.miniset.786O = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

write.csv(res.lfc.DMSO.day0.miniset.786O,"Processed_data/786O_DMSO_vs_Day0_CRISPRa_HIF2a_library.csv")
write.csv(res.lfc.PT.day0.miniset.786O,"Processed_data/786O_PT_vs_Day0_CRISPRa_HIF2a_library.csv")
write.csv(res.lfc.PT.DMSO.miniset.786O,"Processed_data/786O_PT_vs_DMSO_CRISPRa_HIF2a_library.csv")
##################################################################################################################
##################################################################################################################

# Volcano plots

num_genes_highlight_sens <- 15
num_genes_highlight_res <- 15

plot_file <- res.lfc.DMSO.day0.comb 
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(-3,3))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+ 
  ggtitle("DMSO d20 vs d0")+
  guides(color = FALSE)
ggsave(gg, filename = "Figures/FigS1C/Volcano_DMSO_d0_comb_TUHR4TKB.pdf", width = 10, height = 10, dpi = 420)

plot_file <- res.lfc.PT.day0.comb 
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(-3,3))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+  guides(color = FALSE)+
  ggtitle("PT d20 vs d0")
ggsave(gg, filename = "Figures/FigS1C/Volcano_PT_d0_comb_TUHR4TKB.pdf", width = 10, height = 10, dpi = 420)


plot_file <- res.lfc.PT.DMSO.comb 
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(-3,3))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+ guides(color = FALSE) +
  ggtitle("PT d20 vs DMSO d20")
ggsave(gg, filename = "Figures/FigS1C/Volcano_PT_DMSO_comb_TUHR4TKB.pdf", width = 10, height = 10, dpi = 420)

plot_file <- res.lfc.DMSO.day0.comb.OSRC2 
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(-3,3))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+ 
  ggtitle("DMSO d20 vs d0")+
  guides(color = FALSE)
ggsave(gg, filename = "Figures/Fig1C/Volcano_DMSO_d0_comb_OSRC2.pdf", width = 10, height = 10, dpi = 420)

plot_file <- res.lfc.PT.day0.comb.OSRC2
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(-3,3))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+  guides(color = FALSE)+
  ggtitle("PT d20 vs d0")
ggsave(gg, filename = "Figures/Fig1C/Volcano_PT_d0_comb_OSRC2.pdf", width = 10, height = 10, dpi = 420)


plot_file <- res.lfc.PT.DMSO.comb.OSRC2 
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7), limits = c(-3,3))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+ guides(color = FALSE) +
  ggtitle("PT d20 vs DMSO d20")
ggsave(gg, filename = "Figures/Fig1C/Volcano_PT_DMSO_comb_OSRC2.pdf", width = 10, height = 10, dpi = 420)


plot_file <- res.lfc.DMSO.day0.miniset.OSRC2
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11), limits = c(-5,5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+ 
  ggtitle("DMSO d20 vs d0")+
  guides(color = FALSE)
ggsave(gg, filename = "Figures/Fig1G/Volcano_DMSO_d0_HIF2a_subpool_OSRC2.pdf", width = 10, height = 10, dpi = 420)

plot_file <- res.lfc.PT.day0.miniset.OSRC2
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11), limits = c(-5,5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+  guides(color = FALSE)+
  ggtitle("PT d20 vs d0")
ggsave(gg, filename = "Figures/Fig1G/Volcano_PT_d0_HIF2a_subpool_OSRC2.pdf", width = 10, height = 10, dpi = 420)


plot_file <- res.lfc.PT.DMSO.miniset.OSRC2
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(-2,2))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7))+ guides(color = FALSE) +
  ggtitle("PT d20 vs DMSO d20")
ggsave(gg, filename = "Figures/Fig1G/Volcano_PT_DMSO_HIF2a_subpool_OSRC2.pdf", width = 10, height = 10, dpi = 420)


plot_file <- res.lfc.DMSO.day0.miniset.TUHR
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11), limits = c(-5,5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7))+ 
  ggtitle("DMSO d20 vs d0")+
  guides(color = FALSE)
ggsave(gg, filename = "Figures/FigS1G/Volcano_DMSO_d0_HIF2a_subpool_TUHR4TKB.pdf", width = 10, height = 10, dpi = 420)

plot_file <- res.lfc.PT.day0.miniset.TUHR
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 11), limits = c(-5,5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+  guides(color = FALSE)+
  ggtitle("PT d20 vs d0")
ggsave(gg, filename = "Figures/FigS1G/Volcano_PT_d0_HIF2a_subpool_TUHR4TKB.pdf", width = 10, height = 10, dpi = 420)


plot_file <- res.lfc.PT.DMSO.miniset.TUHR
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(-2,2))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0,5))+ guides(color = FALSE) +
  ggtitle("PT d20 vs DMSO d20")
ggsave(gg, filename = "Figures/FigS1G/Volcano_PT_DMSO_HIF2a_subpool_TUHR4TKB.pdf", width = 10, height = 10, dpi = 420)


plot_file <- res.lfc.DMSO.day0.miniset.786O
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(-1,1.5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0,5))+ 
  guides(color = FALSE)+ggtitle("DMSO d20 vs d0")
ggsave(gg, filename = "Figures/FigS3C/Volcano_DMSOd20_vs_d0_HIF2a_subpool_786O.pdf", width = 10, height = 10, dpi = 420)

plot_file <- res.lfc.PT.day0.miniset.786O
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.5,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+ggtitle("PT d20 vs d0")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), limits = c(-1,1.5))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3), limits = c(0,3))+  guides(color = FALSE)
ggsave(gg, filename = "Figures/FigS3C/Volcano_PT_d0_HIF2a_subpool_786O.pdf", width = 10, height = 10, dpi = 420)


plot_file <- res.lfc.PT.DMSO.miniset.786O
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                 slice_head(n=num_genes_highlight_sens))$Gene,"red",
                          ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% 
                                                        filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene)) %>% 
                                                        slice_tail(n=num_genes_highlight_res))$Gene
                                 ,"blue","gray"))
gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene, color = color)) +
  geom_point(size = 6) + theme_classic() + 
  ggrepel::geom_text_repel(data = plot_file %>% filter(color %in% c("red","blue")), 
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           force = 1,
                           segment.alpha = .5,
                           segment.size = 0.8, size = 7, segment.color = 'black') +
  theme(axis.text.x = element_text(size = 27),
        axis.text.y = element_text(size = 27),
        axis.title.x = element_text(size = 27),
        axis.title.y = element_text(size = 27),
        plot.title = element_text(size = 30, face = "bold"),
        legend.position="none")+ggtitle("PT d20 vs DMSO d20")+
  scale_color_manual(values=c("blue","gray","red"))+ xlab("Average fold change(log2)") + ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 9), limits = c(-2,2))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7))+ guides(color = FALSE)
ggsave(gg, filename = "Figures/FigS3C/Volcano_PT_vs_DMSO_HIF2a_subpool_786O.pdf", width = 10, height = 10, dpi = 420)

