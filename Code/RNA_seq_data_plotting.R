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

# Plotting data
## Read data
res.lfc.PT.DMSO.OSRC2.24 <- read.csv("Processed_data/RNAseq_OSRC2_24HR_PT_DMSO.csv") %>%
  na.omit() %>% 
  dplyr::select(Gene=symbol, log2FoldChange,pvalue) %>% 
  as.data.frame()
res.lfc.PT.DMSO.OSRC2.48<- read.csv("Processed_data/RNAseq_OSRC2_48HR_PT_DMSO.csv") %>%
  na.omit() %>% 
  dplyr::select(Gene=symbol, log2FoldChange,pvalue) %>% 
  as.data.frame()
res.lfc.PT.DMSO.TUHR.24 <- read.csv("Processed_data/RNAseq_TUHR4TKB_24HR_PT_DMSO.csv") %>%
  na.omit() %>% 
  dplyr::select(Gene=symbol, log2FoldChange,pvalue) %>% 
  as.data.frame()
res.lfc.PT.DMSO.TUHR.48 <- read.csv("Processed_data/RNAseq_TUHR4TKB_48HR_PT_DMSO.csv") %>%
  na.omit() %>% 
  dplyr::select(Gene=symbol, log2FoldChange,pvalue) %>% 
  as.data.frame()

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

### sgRNA abundance PT vs DMSO d20 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_comb,sgRNA_abundance_DMSOd20_comb) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.DMSO.comb.OSRC2 = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

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
### sgRNA abundance PT vs DMSO d20 combined
lfc_matrix <- full_join(sgRNA_abundance_PTd20_comb,sgRNA_abundance_DMSOd20_comb) %>% select(-Construct.IDs) %>%
  group_by(geneSymbol) %>%summarise_all("mean") %>%
  column_to_rownames(var = "geneSymbol") %>% t()
indicator <- as.numeric(list(1,1,0,0))
res.lfc.PT.DMSO.comb = run_lm_stats_limma(mat = lfc_matrix, vec = indicator)

##################################################################################################################
##################################################################################################################

# Volcano plots

plot_file <- res.lfc.PT.DMSO.OSRC2.24 %>% na.omit() %>% 
  dplyr::select(Gene=symbol, EffectSize=log2FoldChange,p.value=pvalue) %>% 
  as.data.frame() %>%
  mutate(color = ifelse(Gene %in% (plot_file %>% filter(EffectSize < 0)%>% arrange(p.value) %>% 
                                     slice_head(n=15))$Gene,"red",
                        ifelse(Gene %in% (plot_file %>% filter(EffectSize > 0)%>% arrange(p.value)%>%
                                            slice_head(n=15))$Gene
                               ,"blue","gray")))

gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene)) +
  geom_point(aes(color = plot_file$color, size = 6)) + theme_classic() +
  scale_color_manual(values=c("blue","gray","red"))+
  ggrepel::geom_text_repel(data = plot_file%>% filter(color!="gray"),
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.8,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           segment.size = 0.5, size = 7, segment.color = 'black') +
  xlab("Average fold change(log2)") + 
  ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position="none")
ggsave(gg, filename = "Figures/FigS2A/Volcano_24hr_PTvsDMSO_RNAseq_OSRC2_HISAT2_subread.pdf", width = 14, height = 8, dpi = 320)

plot_file <- res.lfc.PT.DMSO.OSRC2.48 %>% na.omit() %>% 
  dplyr::select(Gene, EffectSize=log2FoldChange,p.value=pvalue) %>% 
  as.data.frame() %>%
  mutate(color = ifelse(Gene %in% (plot_file %>% filter(EffectSize < 0)%>% arrange(p.value) %>% 
                                     slice_head(n=15))$Gene,"red",
                        ifelse(Gene %in% (plot_file %>% filter(EffectSize > 0)%>% arrange(p.value)%>%
                                            slice_head(n=15))$Gene
                               ,"blue","gray")))

gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene)) +
  geom_point(aes(color = plot_file$color, size = 6)) + theme_classic() +
  scale_color_manual(values=c("blue","gray","red"))+
  ggrepel::geom_text_repel(data = plot_file%>% filter(color!="gray"),
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.8,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           segment.size = 0.5, size = 7, segment.color = 'black') +
  xlab("Average fold change(log2)") + 
  ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position="none")
ggsave(gg, filename = "Figures/FigS2A/Volcano_48hr_PTvsDMSO_RNAseq_OSRC2_HISAT2_subread.pdf", width = 14, height = 8, dpi = 320)

plot_file <- res.lfc.PT.DMSO.TUHR.24 %>% na.omit() %>% 
  dplyr::select(Gene, EffectSize=log2FoldChange,p.value=pvalue) %>% 
  as.data.frame() %>%
  mutate(color = ifelse(Gene %in% (plot_file %>% filter(EffectSize < 0)%>% arrange(p.value) %>% 
                                     slice_head(n=15))$Gene,"red",
                        ifelse(Gene %in% (plot_file %>% filter(EffectSize > 0)%>% arrange(p.value)%>%
                                            slice_head(n=15))$Gene
                               ,"blue","gray")))

gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene)) +
  geom_point(aes(color = plot_file$color, size = 6)) + theme_classic() +
  scale_color_manual(values=c("blue","gray","red"))+
  ggrepel::geom_text_repel(data = plot_file%>% filter(color!="gray"),
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.8,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           segment.size = 0.5, size = 7, segment.color = 'black') +
  xlab("Average fold change(log2)") + 
  ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position="none")
ggsave(gg, filename = "Figures/FigS2B/Volcano_24hr_PTvsDMSO_RNAseq_TUHR4TKB_HISAT2_subread.pdf", width = 14, height = 8, dpi = 320)

plot_file <- res.lfc.PT.DMSO.TUHR.48 %>% na.omit() %>% 
  dplyr::select(Gene, EffectSize=log2FoldChange,p.value=pvalue) %>% 
  as.data.frame() %>%
  mutate(color = ifelse(Gene %in% (plot_file %>% filter(EffectSize < 0)%>% arrange(p.value) %>% 
                                     slice_head(n=10))$Gene,"red",
                        ifelse(Gene %in% (plot_file %>% filter(EffectSize > 0)%>% arrange(p.value)%>%
                                            slice_head(n=10))$Gene
                               ,"blue","gray")))

gg <- ggplot(data = plot_file, aes(x = EffectSize, y = -log10(p.value), label = Gene)) +
  geom_point(aes(color = plot_file$color, size = 6)) + 
  geom_point(data=plot_file%>% filter(Gene=="CCND1"),aes(size = 6),colour="black",pch=21) +
  theme_classic() +
  scale_color_manual(values=c("blue","gray","red"))+
  ggrepel::geom_text_repel(data = plot_file%>% filter(color!="gray"|Gene=="CCND1"),
                           mapping = aes(x = EffectSize, y = -log10(p.value), 
                                         label = Gene), 
                           min.segment.length = 0.8,
                           max.overlaps = 10000,
                           box.padding = 0.5,
                           segment.size = 0.5, size = 7, segment.color = 'black') +
  xlab("Average fold change(log2)") + 
  ylab("Average p-value(-log10)")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 7))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position="none")
ggsave(gg, filename = "Figures/FigS2B/Volcano_48hr_PTvsDMSO_RNAseq_TUHR4TKB_HISAT2_subread.pdf", width = 14, height = 8, dpi = 320)

##################################################################################################################
##################################################################################################################

# Scatterplot

plot_file <- res.lfc.PT.DMSO.comb %>% inner_join(res.lfc.PT.DMSO.TUHR.48) %>% filter(EffectSize >0 & log2FoldChange <0) 
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% filter(EffectSize >= 0.5 & p.value <= 0.05 & log2FoldChange <= -0.5 & pvalue <= 0.05)%>%
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene))%>%
                                                 filter(!grepl("NO_SITE",Gene)))$Gene,"gray","#A100C2")

gg <- ggplot(data = plot_file, aes(x = EffectSize, y = log2FoldChange, label = Gene)) +
  geom_point(aes(color = plot_file$color), size = 6) + theme_classic() +
  ggrepel::geom_text_repel(data = plot_file%>% filter(EffectSize >= 0.5 & p.value <= 0.05 & log2FoldChange <= -0.5 & pvalue <= 0.05),
                           mapping = aes(x = EffectSize, y = log2FoldChange,
                                         label = Gene),
                           min.segment.length = 0,
                           max.overlaps = 10000,
                           segment.size = 0.5, size = 10, segment.color = 'black') +
  scale_color_manual(values=c("gray","#A100C2"),
                     labels = c("Doesn't meet cutoff","Meets cutoff"))+
  xlab("CRISPRa PT2399 vs DMSO average log2 fold change") +
  ylab("RNAseq 48hr log2 foldChange PT2399 vs DMSO")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position="none")
ggsave(gg, filename = "Figures/FigS1E/Scatterplot_CRISPRa_enriched_PTvsDMSO_48hr_RNAseq_HISAT2_TUHR4TKB.pdf", width = 8, height = 8, dpi = 420)

plot_file <- res.lfc.PT.DMSO.comb.OSRC2 %>% inner_join(res.lfc.PT.DMSO.OSRC2.24) %>% filter(EffectSize >0 & log2FoldChange <0) 
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% filter(EffectSize >= 0.5 & p.value <= 0.05 & log2FoldChange <= -0.5 & pvalue <= 0.05)%>%
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene))%>%
                                                 filter(!grepl("NO_SITE",Gene)))$Gene,"#A100C2","gray")

gg <- ggplot(data = plot_file, aes(x = EffectSize, y = log2FoldChange, label = Gene)) +
  geom_point(aes(color = plot_file$color),size = 6) + theme_classic() +
  ggrepel::geom_text_repel(data = plot_file %>% filter(EffectSize >= 0.5 & p.value <= 0.05 & log2FoldChange <= -0.5 & pvalue <= 0.05),
                           mapping = aes(x = EffectSize, y = log2FoldChange,
                                         label = Gene),
                           min.segment.length = 0,
                           box.padding = 0.7,
                           max.overlaps = 10000,
                           segment.size = 0.8, size = 10, segment.color = 'black') +
  scale_color_manual(values=c("#A100C2","gray"),
                     labels = c("Doesn't meet cutoff","Meets cutoff"))+
  xlab("CRISPRa PT2399 vs DMSO average log2 fold change") +
  ylab("RNAseq 24hr log2 foldChange PT2399 vs DMSO")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
  guides(color = "none") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
ggsave(gg, filename = "Figures/Fig1D/Scatterplot_CRISPRa_PTvsDMSO_24hr_RNAseq_HISAT2_OSRC2.pdf", width = 8, height = 8, dpi = 420)

plot_file <- res.lfc.PT.DMSO.comb %>% inner_join(res.lfc.PT.DMSO.TUHR.24) %>% filter(EffectSize >0 & log2FoldChange <0) 
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% filter(EffectSize >= 0.5 & p.value <= 0.05 & log2FoldChange <= -0.5 & pvalue <= 0.05)%>%
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene))%>%
                                                 filter(!grepl("NO_SITE",Gene)))$Gene,"#A100C2","gray")

gg <- ggplot(data = plot_file, aes(x = EffectSize, y = log2FoldChange, label = Gene)) +
  geom_point(aes(color = plot_file$color), size = 6) + theme_classic() +
  ggrepel::geom_text_repel(data = plot_file%>% filter(EffectSize >= 0.5 & p.value <= 0.05 & log2FoldChange <= -0.5 & pvalue <= 0.05),
                           mapping = aes(x = EffectSize, y = log2FoldChange,
                                         label = Gene),
                           min.segment.length = 0,
                           max.overlaps = 10000,
                           segment.size = 0.8, size = 10, segment.color = 'black') +
  scale_color_manual(name = "Cutoffs(LFC, pvalue) - RNAseq(0.5,0.05),CRISPRa(-0.5,0.05)",
                     values=c("#A100C2","gray"),
                     labels = c("Doesn't meet cutoff","Meets cutoff"))+
  xlab("CRISPRa PT2399 vs DMSO average log2 fold change") +
  ylab("RNAseq 24hr log2 foldChange PT2399 vs DMSO")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position="none")
ggsave(gg, filename = "Figures/Fig1E/Scatterplot_CRISPRa_PTvsDMSO_24hr_RNAseq_HISAT2_TUHR4TKB.pdf", width = 8, height = 8, dpi = 420)


plot_file <- res.lfc.PT.DMSO.comb.OSRC2 %>% inner_join(res.lfc.PT.DMSO.OSRC2.48) %>% filter(EffectSize >0 & log2FoldChange <0) 
plot_file$color <- ifelse(plot_file$Gene %in% (plot_file %>% arrange(EffectSize) %>% filter(EffectSize >= 0.5 & p.value <= 0.05 & log2FoldChange <= -0.5 & pvalue <= 0.05)%>%
                                                 filter(!grepl(paste(genes_not_to_include, collapse="|"),Gene))%>%
                                                 filter(!grepl("NO_SITE",Gene)))$Gene,"#A100C2","gray")

gg <- ggplot(data = plot_file, aes(x = EffectSize, y = log2FoldChange, label = Gene)) +
  geom_point(aes(color = plot_file$color), size = 6) + theme_classic() +
  ggrepel::geom_text_repel(data = plot_file%>% filter(EffectSize >= 0.5 & p.value <= 0.05 & log2FoldChange <= -0.5 & pvalue <= 0.05),
                           mapping = aes(x = EffectSize, y = log2FoldChange,
                                         label = Gene),
                           min.segment.length = 0,
                           box.padding = 0.7,
                           max.overlaps = 10000,
                           segment.size = 0.8, size = 10, segment.color = 'black') +
  scale_color_manual(values=c("#A100C2","gray"),
                     labels = c("Doesn't meet cutoff","Meets cutoff"))+
  xlab("CRISPRa PT2399 vs DMSO average log2 fold change") +
  ylab("RNAseq 48hr log2 foldChange PT2399 vs DMSO")+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.position="none")
ggsave(gg, filename = "Figures/FigS1D/Scatterplot_CRISPRa_PTvsDMSO_48hr_RNAseq_HISAT2_OSRC2.pdf", width = 8, height = 8, dpi = 420)

##################################################################################################################
##################################################################################################################

# Venn for scatterplot

rnaseq_TUHR4TKB_scatterplot_24hr <- res.lfc.PT.DMSO.comb %>% 
  inner_join(res.lfc.PT.DMSO.TUHR.24) %>% select(Gene, EffectSize,p.value,TUHR4TKB_24hr_PTvsDMSO_LFC=log2FoldChange,TUHR4TKB_24hr_PTvsDMSO_pvalue=pvalue)%>%
  filter(EffectSize >=0.5 & TUHR4TKB_24hr_PTvsDMSO_LFC <=-0.5 & p.value <= 0.05 & TUHR4TKB_24hr_PTvsDMSO_pvalue<= 0.05) 
rnaseq_TUHR4TKB_scatterplot_48hr <- res.lfc.PT.DMSO.comb %>% 
  inner_join(res.lfc.PT.DMSO.TUHR.48) %>%  select(Gene, EffectSize,p.value,TUHR4TKB_48hr_PTvsDMSO_LFC=log2FoldChange,TUHR4TKB_48hr_PTvsDMSO_pvalue=pvalue)%>%
  filter(EffectSize >=0.5 & TUHR4TKB_48hr_PTvsDMSO_LFC <=-0.5 & p.value <= 0.05 & TUHR4TKB_48hr_PTvsDMSO_pvalue<= 0.05) 
rnaseq_OSRC2_scatterplot_24hr <- res.lfc.PT.DMSO.comb.OSRC2 %>% 
  inner_join(res.lfc.PT.DMSO.OSRC2.24) %>% select(Gene, EffectSize,p.value,OSRC2_24hr_PTvsDMSO_RNAseq=log2FoldChange,OSRC2_24hr_PTvsDMSO_pvalue=pvalue)%>%
  filter(EffectSize >=0.5 & OSRC2_24hr_PTvsDMSO_RNAseq <=-0.5 & p.value <= 0.05 & OSRC2_24hr_PTvsDMSO_pvalue<= 0.05) 
rnaseq_OSRC2_scatterplot_48hr <- res.lfc.PT.DMSO.comb.OSRC2 %>% 
  inner_join(res.lfc.PT.DMSO.OSRC2.48) %>% select(Gene, EffectSize,p.value,OSRC2_48hr_PTvsDMSO_RNAseq=log2FoldChange,OSRC2_48hr_PTvsDMSO_pvalue=pvalue)%>%
  filter(EffectSize >=0.5 & OSRC2_48hr_PTvsDMSO_RNAseq <=-0.5 & p.value <= 0.05 & OSRC2_48hr_PTvsDMSO_pvalue<= 0.05)


x <- list(
  `OSRC2 scatterplot hits RNAseq 24hr-CRISPRa` = rnaseq_OSRC2_scatterplot_24hr$Gene, 
  `TUHR4TKB scatterplot hits RNAseq 24hr-CRISPRa` = rnaseq_TUHR4TKB_scatterplot_24hr$Gene
)
venn <- Venn(x)
data <- process_data(venn, shape_id = "201")
xs <- ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = count,group = id), 
               data = venn_regionedge(data)) +
  # 2. set edge layer
  geom_path(aes(X, Y, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  # 3. set label layer
  geom_text(aes(X, Y, label = name), 
            data = venn_setlabel(data)) +
  coord_equal() +
  theme_void()+
  #scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  theme(legend.position = "none")+
  scale_fill_gradient(low = "#510080", high = "#FFF0FA")
# Add text for elements
grid.text("CCND1", x = 0.48, y = 0.5, gp = gpar(col = "white"))
grid.text("WT1", x = 0.28, y = 0.55, gp = gpar(col = "white"))
grid.text("APOL1", x = 0.28, y = 0.45, gp = gpar(col = "white"))
grid.text("TP73", x = 0.73, y = 0.6, gp = gpar(col = "black"))
grid.text("BHLHE40", x = 0.75, y = 0.53, gp = gpar(col = "black"))
grid.text("MAP3K21", x = 0.75, y = 0.46, gp = gpar(col = "black"))
grid.text("CDKN2C", x = 0.72, y = 0.39, gp = gpar(col = "black"))
ggsave(xs, filename = "Figures/Fig1F/Scatterplot_CRISPRa_PTvsDMSO_24hr_RNAseq_OSRC2_TUHR.pdf", width = 8, height = 8, dpi = 420)
write.csv(apply(venn_region(data),2,as.character),file = "Processed_data/Scatterplot_CRISPRa_PTvsDMSO_24hr_RNAseq_OSRC2_TUHR.csv")

x <- list(
  `OSRC2 scatterplot hits RNAseq 48hr-CRISPRa` = rnaseq_OSRC2_scatterplot_48hr$Gene, 
  `TUHR4TKB scatterplot hits RNAseq 48hr-CRISPRa` = rnaseq_TUHR4TKB_scatterplot_48hr$Gene
)
venn <- Venn(x)
data <- process_data(venn, shape_id = "201")
xs <- ggplot() +
  # 1. region count layer
  geom_polygon(aes(X, Y, fill = count, group = id), 
               data = venn_regionedge(data)) +
  # 2. set edge layer
  geom_path(aes(X, Y, group = id), 
            data = venn_setedge(data), 
            show.legend = FALSE) +
  # 3. set label layer
  geom_text(aes(X, Y, label = name), 
            data = venn_setlabel(data)) +
  coord_equal() +
  theme_void()+
  #scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  theme(legend.position = "none")+
  scale_fill_gradient(low = "#510080", high = "#FFF0FA")
# Add text for elements
grid.text("CCND1", x = 0.5, y = 0.5, gp = gpar(col = "white"))

grid.text("WT1", x = 0.265, y = 0.47, gp = gpar(col = "white"))
grid.text("APOL1", x = 0.28, y = 0.38, gp = gpar(col = "white"))
grid.text("NR4A1", x = 0.27, y = 0.56, gp = gpar(col = "white"))
grid.text("MAP3K8", x = 0.29, y = 0.65, gp = gpar(col = "white"))

grid.text("MYBL1", x = 0.71, y = 0.35, gp = gpar(col = "black"))
grid.text("MAP3K21", x = 0.72, y = 0.4, gp = gpar(col = "black"))
grid.text("SSC5D", x = 0.725, y = 0.45, gp = gpar(col = "black"))
grid.text("CCNE1", x = 0.73, y = 0.5, gp = gpar(col = "black"))
grid.text("TP73", x = 0.73, y = 0.55, gp = gpar(col = "black"))
grid.text("CDKN2C", x = 0.725, y = 0.6, gp = gpar(col = "black"))
grid.text("BHLHE40", x = 0.72, y = 0.65, gp = gpar(col = "black"))
ggsave(xs, filename = "Figures/FigS1F/Scatterplot_CRISPRa_PTvsDMSO_48hr_RNAseq_HISAT2_OSRC2_TUHR.pdf", width = 8, height = 8, dpi = 420)
write.csv(apply(venn_region(data),2,as.character),file = "Processed_data/Scatterplot_CRISPRa_PTvsDMSO_48hr_RNAseq_HISAT2_OSRC2_TUHR.csv")

