rm(list=ls())
library("DESeq2")

setwd(file.path("//salt","wheat_rnaseq", "PeterBuchner_RNA", "Watkins_lines_root_exp", "Analysis", "HISAT_FeatureCounts_DESeq", "5varieties"))
counts <- read.table("counts.tab", row.names = 1, sep="\t", comment="", as.is=T)

colnames(counts) <- c("W01_highN_A","W01_highN_B","W01_highN_C","W01_lowN_A", "W01_lowN_B", "W01_lowN_C",
                      "W02_highN_A","W02_highN_B","W02_highN_C","W02_lowN_A", "W02_lowN_B", "W02_lowN_C",
                      "W05_highN_A","W05_highN_B","W05_highN_C","W05_lowN_A", "W05_lowN_B", "W05_lowN_C",
                      "W10_highN_A","W10_highN_B","W10_highN_C","W10_lowN_A", "W10_lowN_B", "W10_lowN_C",
                      "Paragon_highN_A","Paragon_highN_B","Paragon_highN_C","Paragon_lowN_A", "Paragon_lowN_B", "Paragon_lowN_C")

sampleInfo <- data.frame(row.names= colnames(counts),
                         variety=factor(c(rep("W01", 6), rep("W02", 6), rep("W05", 6), rep("W10", 6), rep("Paragon", 6)), levels=c("Paragon", "W01", "W02", "W05", "W10")),
                         condition=factor(rep(c(rep("highN",3),rep("lowN",3)),5), levels=c("lowN", "highN")))

#Separate highN and lowN and drop Paragon to identify genes with differential expression in the lines with contrasting phenotypes
library("tibble")
library("dplyr")
highN_watkins_counts <- select(counts, contains("highN")) %>% select(-starts_with("Paragon"))
lowN_watkins_counts <- select(counts, contains("lowN")) %>% select(-starts_with("Paragon"))

highNInfo <- rownames_to_column(sampleInfo) %>% filter(condition=="highN" & variety != "Paragon") %>% column_to_rownames() %>% select(variety)
highNInfo$variety <- relevel(highNInfo$variety, "W01") %>% droplevels()
lowNInfo <- rownames_to_column(sampleInfo) %>% filter(condition=="lowN" & variety != "Paragon") %>% column_to_rownames() %>% select(variety)
lowNInfo$variety <- relevel(lowNInfo$variety, "W01") %>% droplevels()

#ddsMat <- DESeqDataSetFromMatrix(countData = counts, sampleInfo, ~ variety + condition)
sampleInfo$group <- factor(paste(sampleInfo$variety, sampleInfo$condition, sep = "_"))
ddsMat <- DESeqDataSetFromMatrix(countData = counts, sampleInfo, ~ 0 + group)
dds <- ddsMat[rowSums(counts(ddsMat)) > 1,]
rld <- rlog(dds)

#Repeat for highN and lowN
ddsMat_highN <- DESeqDataSetFromMatrix(countData = highN_watkins_counts, highNInfo, ~ 0 + variety)
dds_highN <- ddsMat_highN[rowSums(counts(ddsMat_highN)) > 1,]

ddsMat_lowN <- DESeqDataSetFromMatrix(countData = lowN_watkins_counts, lowNInfo, ~ 0 + variety)
dds_lowN <- ddsMat_lowN[rowSums(counts(ddsMat_lowN)) > 1,]

#Use DESeq and plotly to create an interactive PCA 
library("plotly")
#PCAplot <- plotPCA(rld, intgroup = c("variety","condition"))
PCAplot <- plotPCA(rld, intgroup = "group")
p <- ggplotly(PCAplot)
p <- style(p, hoverinfo="name")
p

#####################################################################################
#Genes with different means in highN vs lowN
#Paragon included in the analysis
#####################################################################################

dds <- DESeq(dds)
resultsNames(dds)

res_N <- results(dds, contrast = c(1/5,-1/5,1/5,-1/5,1/5,-1/5,1/5,-1/5,1/5,-1/5), alpha=0.05)
#rownames(res_N) <- rownames(res_N) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_N)
degs_N <- subset(res_N, padj < 0.05)
dim(degs_N)
degs_N <- degs_N[order(degs_N$padj),]
#write.table(degs, file="paragon_degs.tab", quote = F, sep = "\t", row.names = T)
#write(rownames(degs_N), "degs_N_ids.txt")

#Too general, filter applied to genes that are differentially expressed between highN and lowN in all varieties
res_N_W01 <- results(dds, contrast = list("groupW01_highN", "groupW01_lowN"), alpha=0.05)
#rownames(res_N_W01) <- rownames(res_N_W01) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_N_W01)
degs_N_W01 <- subset(res_N_W01, padj < 0.05)
dim(degs_N_W01)
degs_N_W01 <- degs_N_W01[order(degs_N_W01$padj),]
#write(rownames(degs_N_W01), "degs_N_W01_ids.txt")
#write.xlsx(degs_N_W01,
#           "HIvLO_N_degs.xlsx",
#           sheetName = "W01_DEGs",
#           row.names = T,
#           append = T)

res_N_W02 <- results(dds, contrast = list("groupW02_highN", "groupW02_lowN"), alpha=0.05)
#rownames(res_N_W02) <- rownames(res_N_W02) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_N_W02)
degs_N_W02 <- subset(res_N_W02, padj < 0.05)
dim(degs_N_W02)
degs_N_W02 <- degs_N_W02[order(degs_N_W02$padj),]
#write(rownames(degs_N_W02), "degs_N_W02_ids.txt")
#write.xlsx(degs_N_W02,
#           "HIvLO_N_degs.xlsx",
#           sheetName = "W02_DEGs",
#           row.names = T,
#           append = T)

res_N_Paragon <- results(dds, contrast = list("groupParagon_highN", "groupParagon_lowN"), alpha=0.05)
#rownames(res_N_Paragon) <- rownames(res_N_Paragon) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_N_Paragon)
degs_N_Paragon <- subset(res_N_Paragon, padj < 0.05)
dim(degs_N_Paragon)
degs_N_Paragon <- degs_N_Paragon[order(degs_N_Paragon$padj),]
#write(rownames(degs_N_Paragon), "degs_N_Paragon_ids.txt")
write.xlsx(degs_N_Paragon,
           "HIvLO_N_degs.xlsx",
           sheetName = "Paragon_DEGs",
           row.names = T,
           append = T)

res_N_W05 <- results(dds, contrast = list("groupW05_highN", "groupW05_lowN"), alpha=0.05)
#rownames(res_N_W05) <- rownames(res_N_W05) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_N_W05)
degs_N_W05 <- subset(res_N_W05, padj < 0.05)
dim(degs_N_W05)
degs_N_W05 <- degs_N_W05[order(degs_N_W05$padj),]
#write(rownames(degs_N_W05), "degs_N_W05_ids.txt")
#write.xlsx(degs_N_W05,
#           "HIvLO_N_degs.xlsx",
#           sheetName = "W05_DEGs",
#           row.names = T,
#           append = T)

res_N_W10 <- results(dds, contrast = list("groupW10_highN", "groupW10_lowN"), alpha=0.05)
#rownames(res_N_W10) <- rownames(res_N_W10) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_N_W10)
degs_N_W10 <- subset(res_N_W10, padj < 0.05)
dim(degs_N_W10)
degs_N_W10 <- degs_N_W10[order(degs_N_W10$padj),]
#write(rownames(degs_N_W10), "degs_N_W10_ids.txt")
#write.xlsx(degs_N_W10,
#           "HIvLO_N_degs.xlsx",
#           sheetName = "W10_DEGs",
#           row.names = T,
#           append = T)

#Identify genes with N effect that are differentially expressed between highN and lowN in all varieties
library("magrittr")
filtered_degs_N_watkins <- rownames(degs_N_W01) %>%
                            intersect(rownames(degs_N_W02)) %>%
                            intersect(rownames(degs_N_W05)) %>%
                            intersect(rownames(degs_N_W10)) %>% 
                            intersect(rownames(degs_N))

filtered_degs_N_all <- rownames(degs_N_Paragon) %>%
                        intersect(filtered_degs_N_watkins)

filtered_degs_N_watkins_noParagon <- setdiff(filtered_degs_N_watkins, filtered_degs_N_all)

# filtered_degs_N_down <- rownames(degs_N_Paragon[degs_N_Paragon$log2FoldChange<0,]) %>%
#                         intersect(rownames(degs_N_W01[degs_N_W01$log2FoldChange<0,])) %>%
#                         intersect(rownames(degs_N_W02[degs_N_W02$log2FoldChange<0,])) %>%
#                         intersect(rownames(degs_N_W05[degs_N_W05$log2FoldChange<0,])) %>%
#                         intersect(rownames(degs_N_W10[degs_N_W10$log2FoldChange<0,])) %>% 
#                         intersect(rownames(degs_N[degs_N$log2FoldChange<0,]))
# 
# filtered_degs_N_up <- rownames(degs_N_Paragon[degs_N_Paragon$log2FoldChange>0,]) %>%
#                       intersect(rownames(degs_N_W01[degs_N_W01$log2FoldChange>0,])) %>%
#                       intersect(rownames(degs_N_W02[degs_N_W02$log2FoldChange>0,])) %>%
#                       intersect(rownames(degs_N_W05[degs_N_W05$log2FoldChange>0,])) %>%
#                       intersect(rownames(degs_N_W10[degs_N_W10$log2FoldChange>0,])) %>% 
#                       intersect(rownames(degs_N[degs_N$log2FoldChange>0,]))

#Identify genes with N effect only differentially expressed in W01/W02 or W05/W10
filtered_degs_N_W01W02 <- rownames(degs_N_W01) %>%
                          intersect(rownames(degs_N_W02)) %>%
                          intersect(rownames(degs_N))

filtered_degs_N_W01W02_noW05W10 <- filtered_degs_N_W01W02 %>%
                                    setdiff(rownames(degs_N_W05)) %>%
                                    setdiff(rownames(degs_N_W10))
                          
filtered_degs_N_W05W10 <- rownames(degs_N_W05) %>%
                          intersect(rownames(degs_N_W10)) %>%
                          intersect(rownames(degs_N))

filtered_degs_N_W05W10_noW01W02 <- filtered_degs_N_W05W10 %>% 
                                    setdiff(rownames(degs_N_W01)) %>%
                                    setdiff(rownames(degs_N_W02))

#Final set of filtered DEGs used for the homeolog analysis will be the union of those genes that are:
# DE in all varieties (W01, W02, W05, W10 and Paragon)
# DE in W01 and W02 but NOT in W05 and W10
# DE in W05 and W10 but NOT in W01 and W02
filtered_degs_N <- filtered_degs_N_all %>%
                    union(filtered_degs_N_W01W02_noW05W10) %>%
                    union(filtered_degs_N_W05W10_noW01W02)

#Identify genes only DE in Paragon
filtered_degs_N_ParagonOnly <- rownames(degs_N_Paragon) %>%
  intersect(rownames(degs_N)) %>% 
  setdiff(rownames(degs_N_W01)) %>%
  setdiff(rownames(degs_N_W02)) %>%
  setdiff(rownames(degs_N_W05)) %>%
  setdiff(rownames(degs_N_W10))
  

#Save DEGs genes in separate sheets depending on their pattern
options(java.parameters = "- Xmx1024m")
library(xlsx)
#Save all the DEGs between high and low N in the main sheet
write.xlsx(degs_N,
           "HIvLO_N_degs.xlsx",
           sheetName = "mean_across_varieties_DE",
           row.names = T,
           append = F)

#Save genes DE in all varieties
write.xlsx(cbind(degs_N_W01[filtered_degs_N_all,],degs_N_W02[filtered_degs_N_all,],
                 degs_N_W05[filtered_degs_N_all,], degs_N_W10[filtered_degs_N_all,],
                 degs_N_Paragon[filtered_degs_N_all,]),
           "HIvLO_N_degs.xlsx",
           sheetName = "all_varieties_DE",
           row.names = T,
           append = T)

#Save genes DE in watkins varieties
filtered_degs_N_watkins_only <- setdiff(filtered_degs_N_watkins, filtered_degs_N_all)
write.xlsx(cbind(degs_N_W01[filtered_degs_N_watkins_only,],degs_N_W02[filtered_degs_N_watkins_only,],
                 degs_N_W05[filtered_degs_N_watkins_only,], degs_N_W10[filtered_degs_N_watkins_only,],
                 res_N_Paragon[filtered_degs_N_watkins_only,]),
           "HIvLO_N_degs.xlsx",
           sheetName = "watkins_varieties_only_DE",
           row.names = T,
           append = T)

#Save genes DE in W01 and W02 but not in W05 and W10
write.xlsx(cbind(degs_N_W01[filtered_degs_N_W01W02_noW05W10,],degs_N_W02[filtered_degs_N_W01W02_noW05W10,],
                 res_N_W05[filtered_degs_N_W01W02_noW05W10,], res_N_W10[filtered_degs_N_W01W02_noW05W10,]),
           "HIvLO_N_degs.xlsx",
           sheetName = "W01W02DE_noW05W10",
           row.names = T,
           append = T)

#Save genes DE in W05 and W10 but not in W01 and W02
write.xlsx(cbind(degs_N_W05[filtered_degs_N_W05W10_noW01W02,],degs_N_W10[filtered_degs_N_W05W10_noW01W02,],
                 res_N_W01[filtered_degs_N_W05W10_noW01W02,],res_N_W02[filtered_degs_N_W05W10_noW01W02,]),
           "HIvLO_N_degs.xlsx",
           sheetName = "W05W10DE_noW01W02",
           row.names = T,
           append = T)

#Find genes with different fold changes in W01/W02 and W05/W10
filtered_degs_N_all_FC <- cbind(degs_N_W01[filtered_degs_N_all,] %>% as.data.frame() %>% select(log2FoldChange),
      degs_N_W02[filtered_degs_N_all,] %>% as.data.frame() %>% select(log2FoldChange),
      degs_N_W05[filtered_degs_N_all,]%>% as.data.frame() %>% select(log2FoldChange),
      degs_N_W10[filtered_degs_N_all,]%>% as.data.frame() %>% select(log2FoldChange))

filtered_degs_N_all_sameFC <- apply(filtered_degs_N_all_FC, 1,function(gene) gene %>% sign() %>% mean %>% abs() == 1)

#Save genes DE in Paragon only
write.xlsx(cbind(degs_N_Paragon[filtered_degs_N_ParagonOnly,],
                 res_N_W01[filtered_degs_N_ParagonOnly,],res_N_W02[filtered_degs_N_ParagonOnly,],
                 res_N_W05[filtered_degs_N_ParagonOnly,],res_N_W10[filtered_degs_N_ParagonOnly,]),
           "HIvLO_N_degs.xlsx",
           sheetName = "Paragon_only_DE",
           row.names = T,
           append = T)

###Filter based on differential expression of homeologs
source(file.path("//salt","wheat_rnaseq", "NRgeneV1", "ensembl_homeologs_utils.R"))
#homeologs <- parse_iwgsc_homeologs()

#1 Genes with all the homeologs being DE
all_homeologs_degs_N <- filtered_degs_N[sapply(filtered_degs_N,
                                                       function(deg) all_homeologs_appear(deg, filtered_degs_N,F)) %>%
                                                  unlist()]

#1.1 Direction of fold-change (FC of same sign)
all_homeologs_degs_sameFC_N<- sapply(all_homeologs_degs_N, function(deg) genes_same_FC_sign(deg,  res_N))

#Save genes whose homeologs are DE in the same direction
write.xlsx(homeolog_results_as_data_frame(all_homeologs_degs_N[all_homeologs_degs_sameFC_N], res_N),
           "HIvLO_N_filtered_degs.xlsx",
           sheetName = "all_homeologs_deg_sameFC",
           row.names = F,
           append = F)

#Save genes whose homeologs are DE in different directions
#write.xlsx(homeolog_results_as_data_frame(all_homeologs_degs_N[!all_homeologs_degs_sameFC_N], res_N),
#           "filtered_results.xlsx",
#           sheetName = "HIvLO_N_all_homeologs_deg_diffFC",
#           row.names = F,
#           append = T)

#DEGs with no homeologs
singleton_degs_N <- filtered_degs_N[sapply(filtered_degs_N,
                                                   function(deg) is_singleton(deg))]

write.xlsx(degs_N[singleton_degs_N, c("baseMean", "log2FoldChange", "padj")],
           "HIvLO_N_filtered_degs.xlsx",
           sheetName = "singleton_degs",
           row.names = T,
           append = T)

#DEGs with more than one homeologs (but NOT all) DE
#These genes are divided categories based on the expression of their homeologs:
# 1. None of the non-DE homeologs are expressed (were filtered out during processing)
# 2. The baseMean difference between the DEG with the lowest expression and the non-DE homeologs is more than one order of magnitude (>10x)
any_homeologs_degs_N <- filtered_degs_N[sapply(filtered_degs_N,
                                                       function(deg) any_homeologs_appear(deg, filtered_degs_N,F)) %>%
                                                  unlist()]
multiple_homeologs_degs_N <- setdiff(any_homeologs_degs_N, all_homeologs_degs_N)

#Save all the DEGs with multiple homeologs DE
write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs_N, res_N),
           "HIvLO_N_filtered_degs.xlsx",
           sheetName = "multiple_homeologs_deg",
           row.names = F,
           append = T)

#Check the fold-change of the DE homeologs
multiple_homeologs_degs_sameFC_N<- sapply(multiple_homeologs_degs_N, function(deg) {
                                                                        gene_homeologs <- find_homeologs(deg)
                                                                        deg_homeologs <- c(deg, gene_homeologs[gene_homeologs %in% any_homeologs_degs_N])
                                                                        genes_same_FC_sign(deg_homeologs,  res_N)})

multiple_homeologs_degs_nondegs_low_N <- sapply(multiple_homeologs_degs_N,
                                                    function(gene){
                                                      gene_homeologs <- find_homeologs(gene)
                                                      #deg_homeologs <- gene_homeologs[gene_homeologs %in% multiple_homeologs_degs_N]
                                                      deg_homeologs <- c(gene, gene_homeologs[gene_homeologs %in% any_homeologs_degs_N])
                                                      homeolog_results <- retrieve_gene_results(c(gene,gene_homeologs), res_N)
                                                      min_deg_baseMean <- homeolog_results[deg_homeologs, "baseMean"] %>% min()
                                                      nondeg_results <- homeolog_results[! rownames(homeolog_results) %in% deg_homeologs,]
                                                      FC2DEG <- min_deg_baseMean/nondeg_results[, "baseMean"]
                                                      #print(gene)
                                                      if(FC2DEG %>% is.na() %>% sum()==length(c(gene,gene_homeologs))-length(deg_homeologs) || all(FC2DEG > 10, na.rm = T) ) return(T) else return(F)})
multiple_homeologs_degs_nondegs_low_N <- multiple_homeologs_degs_N[multiple_homeologs_degs_nondegs_low_N]

#Save DEGs with multiple homeologs DE kept after selecting the lowly expressed non-DE homeologs
write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs_nondegs_low_N, res_N),
           "HIvLO_N_filtered_degs.xlsx",
           sheetName = "multiple_homeo_deg_low_nondegs",
           row.names = F,
           append = T)

#DEGs with no DE homeologs but whose expression is much higher than that of the homeologs
#Genes will be kept if:
# 1. None of the other homeologs are expressed (were filtered out during processing)
# 2. The baseMean difference between the DEG and the homeologs is more than one order of magnitude (>10x)
#no_homeologs_deg_N <- sapply(filtered_degs_N,
#                                 function(deg) any_homeologs_appear(deg, filtered_degs_N,T)) %>% unlist()
no_homeologs_deg_N <- filtered_degs_N[!filtered_degs_N %in% c(any_homeologs_degs_N, singleton_degs_N)]

#Save all the degs with no homeologs DE
write.xlsx(homeolog_results_as_data_frame(no_homeologs_deg_N, res_N),
           "HIvLO_N_filtered_degs.xlsx",
           sheetName = "one_homeolog_deg",
           row.names = F,
           append = T)

one_deg_low_homeologs_N <- sapply(no_homeologs_deg_N,
                                      function(gene){
                                        gene_homeologs <- find_homeologs(gene)
                                        homeolog_results <- retrieve_gene_results(c(gene,gene_homeologs), res_N) 
                                        FC2DEG <- homeolog_results[gene, "baseMean"]/homeolog_results[, "baseMean"]
                                        if(FC2DEG %>% is.na() %>% sum()==length(FC2DEG)-1 || all(FC2DEG > 10, na.rm = T) ) return(T) else return(F)})
one_deg_low_homeologs_N <- no_homeologs_deg_N[one_deg_low_homeologs_N]

write.xlsx(homeolog_results_as_data_frame(one_deg_low_homeologs_N, res_N),
           "HIvLO_N_filtered_degs.xlsx",
           sheetName = "one_homeolog_deg_low_nondegs",
           row.names = F,
           append = T)


#####################################################################################
#Genes with different means in W01/W02 vs W05/W10
#highN and lowN analysed separately
#####################################################################################

#####################################################################################
#####################################################################################
#highN
dds_highN <- DESeq(dds_highN)
resultsNames(dds_highN)

res_highN <- results(dds_highN, contrast = c(-0.5,-0.5,0.5,0.5), alpha=0.05)
#rownames(res_highN) <- rownames(res_highN) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_highN)
degs_highN <- subset(res_highN, padj < 0.05)
dim(degs_highN)
degs_highN <- degs_highN[order(degs_highN$padj),]

#degs_highN_up <- rownames(degs_highN[degs_highN$log2FoldChange>0,])
#degs_highN_down <- rownames(degs_highN[degs_highN$log2FoldChange<0,])


#More stringent condition to filter long list of DEGs. Gene will also need to be...
# Differentially expressed between W05 and W01
# and differentially expressed between W05 and W02
# and differentially expressed between W10 and W01
# and differentially expressed between W10 and W02

res_W05vsW01_highN <- results(dds_highN, contrast = c(-1,0,1,0), alpha=0.05)
#rownames(res_W05vsW01_highN) <- rownames(res_W05vsW01_highN) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_W05vsW01_highN)
degs_W05vsW01_highN <- subset(res_W05vsW01_highN, padj < 0.05)
dim(degs_W05vsW01_highN)
degs_W05vsW01_highN <- degs_W05vsW01_highN[order(degs_W05vsW01_highN$padj),]

res_W05vsW02_highN <- results(dds_highN, contrast = c(0,-1,1,0), alpha=0.05)
#rownames(res_W05vsW02_highN) <- rownames(res_W05vsW02_highN) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_W05vsW02_highN)
degs_W05vsW02_highN <- subset(res_W05vsW02_highN, padj < 0.05)
dim(degs_W05vsW02_highN)
degs_W05vsW02_highN <- degs_W05vsW02_highN[order(degs_W05vsW02_highN$padj),]

res_W10vsW01_highN <- results(dds_highN, contrast = c(-1,0,0,1), alpha=0.05)
#rownames(res_W10vsW01_highN) <- rownames(res_W10vsW01_highN) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_W10vsW01_highN)
degs_W10vsW01_highN <- subset(res_W10vsW01_highN, padj < 0.05)
dim(degs_W05vsW01_highN)
degs_W10vsW01_highN <- degs_W10vsW01_highN[order(degs_W10vsW01_highN$padj),]

res_W10vsW02_highN <- results(dds_highN, contrast = c(0,-1,0,1), alpha=0.05)
#rownames(res_W10vsW02_highN) <- rownames(res_W10vsW02_highN) %>% gsub("02G", "01G", .)#Change gene ids to official v1.0 identifiers
summary(res_W10vsW02_highN)
degs_W10vsW02_highN <- subset(res_W10vsW02_highN, padj < 0.05)
dim(degs_W10vsW02_highN)
degs_W10vsW02_highN <- degs_W10vsW02_highN[order(degs_W10vsW02_highN$padj),]

# res_W05vsW10_highN <- results(dds_highN, contrast = c(0,0,1,-1), alpha=0.05)
# summary(res_W05vsW10_highN)
# degs_W05vsW10_highN <- subset(res_W05vsW10_highN, padj < 0.05)
# dim(degs_W05vsW10_highN)
# degs_W05vsW10_highN <- degs_W05vsW10_highN[order(degs_W05vsW10_highN$padj),]
# 
# res_W01vsW02_highN <- results(dds_highN, contrast = c(1,-1,0,0), alpha=0.05)
# summary(res_W01vsW02_highN)
# degs_W01vsW02_highN <- subset(res_W01vsW02_highN, padj < 0.05)
# dim(degs_W01vsW02_highN)
# degs_W01vsW02_highN <- degs_W01vsW02_highN[order(degs_W01vsW02_highN$padj),]

filtered_degs_highN <- rownames(degs_highN) %>%
                        intersect(rownames(degs_W05vsW01_highN)) %>%
                        #filtered_degs_highN <- rownames(degs_W05vsW01_highN) %>%
                        intersect(rownames(degs_W05vsW02_highN)) %>% 
                        intersect(rownames(degs_W10vsW01_highN)) %>%
                        intersect(rownames(degs_W10vsW02_highN)) #%>% 
#                        setdiff(rownames(degs_W01vsW02_highN)) %>% 
#                        setdiff(rownames(degs_W05vsW10_highN))

#filtered_degs_highN_up <- rownames(degs_highN[filtered_degs_highN,][degs_highN[filtered_degs_highN,]$log2FoldChange>0,])
#filtered_degs_highN_down <- rownames(degs_highN[filtered_degs_highN,][degs_highN[filtered_degs_highN,]$log2FoldChange<0,])

#Save DEGs genes in separate sheets depending on their pattern
options(java.parameters = "- Xmx1024m")
library(xlsx)

#Save all the DEGs between high and low N in the main sheet
write.xlsx(degs_highN,
           "highN_W01W02vW05W10_degs.xlsx",
           sheetName = "means_DE",
           row.names = T,
           append = F)

#Save genes DE in all varieties
write.xlsx(cbind(degs_W05vsW01_highN[filtered_degs_highN,], degs_W10vsW01_highN[filtered_degs_highN,],
                 degs_W05vsW02_highN[filtered_degs_highN,], degs_W10vsW02_highN[filtered_degs_highN,]),
           "highN_W01W02vW05W10_degs.xlsx",
           sheetName = "varieties_pairwise_DE",
           row.names = T,
           append = T)



###Filter based on differential expression of homeologs
source(file.path("//salt","wheat_rnaseq", "NRgeneV1", "ensembl_homeologs_utils.R"))
#homeologs <- parse_iwgsc_homeologs()

#1 Genes with all the homeologs being DE
all_homeologs_degs_highN <- filtered_degs_highN[sapply(filtered_degs_highN,
                                                 function(deg) all_homeologs_appear(deg, filtered_degs_highN,F)) %>%
                                            unlist()]

#1.1 Direction of fold-change (FC of same sign)
all_homeologs_degs_sameFC_highN<- sapply(all_homeologs_degs_highN, function(deg){
  gene_homeologs <- find_homeologs(deg)
  genes_same_FC_sign(gene_homeologs,  res_highN)})
 


#Save genes whose homeologs are DE in the same direction
#write.xlsx(homeolog_results_as_data_frame(all_homeologs_degs_highN[all_homeologs_degs_sameFC_highN], homeologs, res_highN),
#           "highN_W01W02vW05W10_filtered_degs.xlsx",
#           sheetName = "all_homeologs_deg_sameFC",
#           row.names = F)

#Save genes whose homeologs are DE in different directions
#write.xlsx(homeolog_results_as_data_frame(all_homeologs_degs_highN[!all_homeologs_degs_sameFC_highN], homeologs, res_highN),
#           "filtered_results.xlsx",
#           sheetName = "highN_all_homeologs_deg_diffFC",
#           row.names = F,
#           append = T)

#DEGs with no homeologs
singleton_degs_highN <- filtered_degs_highN[sapply(filtered_degs_highN,
                                                      function(deg) is_singleton(deg))]

#write.xlsx(degs_highN[singleton_degs_highN, c("baseMean", "log2FoldChange", "padj")],
#           "filtered_results.xlsx",
#           sheetName = "highN_singleton_degs",
#           row.names = T,
#           append = T)

#DEGs with more than one homeologs (but NOT all) DE
#Genes will be kept if:
# 1. None of the non-DE homeologs are expressed (were filtered out during processing)
# 2. The baseMean difference between the DEG with the lowest expression and the non-DE homeologs is more than one order of magnitude (>10x)
any_homeologs_degs_highN <- filtered_degs_highN[sapply(filtered_degs_highN,
                                                       function(deg) any_homeologs_appear(deg, filtered_degs_highN,F)) %>%
                                                  unlist()]
multiple_homeologs_degs_highN <- setdiff(any_homeologs_degs_highN, all_homeologs_degs_highN)

#Save all the DEGs with multiple homeologs DE
#write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs_N, res_N),
#           "HIvLO_N_filtered_degs.xlsx",
#           sheetName = "multiple_homeologs_deg",
#           row.names = F,
#           append = T)

#Check the fold-change of the DE homeologs
multiple_homeologs_degs_sameFC_highN<- sapply(multiple_homeologs_degs_highN, function(deg) {
                                                                              gene_homeologs <- find_homeologs(deg)
                                                                              deg_homeologs <- c(deg, gene_homeologs[gene_homeologs %in% any_homeologs_degs_highN])
                                                                              genes_same_FC_sign(deg_homeologs,  res_highN)})

#write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs_highN[multiple_homeologs_degs_sameFC_highN], res_highN),
#           "filtered_results.xlsx",
#           sheetName = "highN_multiple_homeologs_deg",
#           row.names = F,
#           append = T)
#write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs_highN[!multiple_homeologs_degs_sameFC_highN], res_highN),
#           "filtered_results.xlsx",
#           sheetName = "highN_multiple_homeologs_deg",
#           row.names = F,
#           append = T)

multiple_homeologs_degs_nondegs_low_highN <- sapply(multiple_homeologs_degs_highN,
                                                    function(gene){
                                                      gene_homeologs <- find_homeologs(gene)
                                                      deg_homeologs <- c(gene, gene_homeologs[gene_homeologs %in% any_homeologs_degs_highN])
                                                      homeolog_results <- retrieve_gene_results(c(gene,gene_homeologs), res_highN)
                                                      min_deg_baseMean <- homeolog_results[deg_homeologs, "baseMean"] %>% min()
                                                      nondeg_results <- homeolog_results[! rownames(homeolog_results) %in% deg_homeologs,]
                                                      FC2DEG <- min_deg_baseMean/nondeg_results[, "baseMean"]
                                                      #print(gene)
                                                      if(FC2DEG %>% is.na() %>% sum()==length(c(gene,gene_homeologs))-length(deg_homeologs) || all(FC2DEG > 10, na.rm = T) ) return(T) else return(F)})

multiple_homeologs_degs_nondegs_low_highN <- multiple_homeologs_degs_highN[multiple_homeologs_degs_nondegs_low_highN]

#write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs_nondegs_low_highN, homeologs, res_highN),
#           "filtered_results.xlsx",
#           sheetName = "highN_multiple_homeologs_deg",
#           row.names = F,
#           append = T)


#DEGs with no DE homeologs but whose expression is much higher than that of the homeologs
#Genes will be kept if:
# 1. None of the other homeologs are expressed (were filtered out during processing)
# 2. The baseMean difference between the DEG and the homeologs is more than one order of magnitude (>10x)
#no_homeologs_deg_highN <- sapply(filtered_degs_highN,
#                                 function(deg) any_homeologs_appear(deg, homeologs,filtered_degs_highN,T)) %>% unlist()
no_homeologs_deg_highN <- filtered_degs_highN[!filtered_degs_highN %in% c(any_homeologs_degs_highN, singleton_degs_highN)]

#Save all the degs with no homeologs DE
#write.xlsx(homeolog_results_as_data_frame(no_homeologs_deg_N, res_N),
#           "HIvLO_N_filtered_degs.xlsx",
#           sheetName = "one_homeolog_deg",
#           row.names = F,
#           append = T)

one_deg_low_homeologs_highN <- sapply(no_homeologs_deg_highN,
                            function(gene){
                              gene_homeologs <- find_homeologs(gene)
                              homeolog_results <- retrieve_gene_results(c(gene,gene_homeologs), res_highN) 
                              FC2DEG <- homeolog_results[gene, "baseMean"]/homeolog_results[, "baseMean"]
                              if(FC2DEG %>% is.na() %>% sum()==length(FC2DEG)-1 || all(FC2DEG > 10, na.rm = T) ) return(T) else return(F)})
one_deg_low_homeologs_highN <- no_homeologs_deg_highN[one_deg_low_homeologs_highN]

write.xlsx(homeolog_results_as_data_frame(one_deg_low_homeologs_highN, homeologs, res_highN),
           "filtered_results.xlsx",
           sheetName = "highN_one_deg_low_homeologs",
           row.names = F,
           append = T)


#####################################################################################
#####################################################################################
#Repeat for lowN
dds_lowN <- DESeq(dds_lowN)
resultsNames(dds_lowN)

res_lowN <- results(dds_lowN, contrast = c(-0.5,-0.5,0.5,0.5), alpha=0.05)
summary(res_lowN)
degs_lowN <- subset(res_lowN, padj < 0.05)
dim(degs_lowN)
degs_lowN <- degs_lowN[order(degs_lowN$padj),]

#degs_lowN_up <- rownames(degs_lowN[degs_lowN$log2FoldChange>0,])
#degs_lowN_down <- rownames(degs_lowN[degs_lowN$log2FoldChange<0,])

#Filter based on genes being DE both in highN and lowN
#deg_bothN_up <- intersect(degs_highN_up, degs_lowN_up)
#deg_bothN_down <- intersect(degs_highN_down, degs_lowN_down)

res_W05vsW01_lowN <- results(dds_lowN, contrast = c(-1,0,1,0), alpha=0.05)
summary(res_W05vsW01_lowN)
degs_W05vsW01_lowN <- subset(res_W05vsW01_lowN, padj < 0.05)
dim(degs_W05vsW01_lowN)
degs_W05vsW01_lowN <- degs_W05vsW01_lowN[order(degs_W05vsW01_lowN$padj),]

res_W05vsW02_lowN <- results(dds_lowN, contrast = c(0,-1,1,0), alpha=0.05)
summary(res_W05vsW02_lowN)
degs_W05vsW02_lowN <- subset(res_W05vsW02_lowN, padj < 0.05)
dim(degs_W05vsW02_lowN)
degs_W05vsW02_lowN <- degs_W05vsW02_lowN[order(degs_W05vsW02_lowN$padj),]

res_W10vsW01_lowN <- results(dds_lowN, contrast = c(-1,0,0,1), alpha=0.05)
summary(res_W10vsW01_lowN)
degs_W10vsW01_lowN <- subset(res_W10vsW01_lowN, padj < 0.05)
dim(degs_W05vsW01_lowN)
degs_W10vsW01_lowN <- degs_W10vsW01_lowN[order(degs_W10vsW01_lowN$padj),]

res_W10vsW02_lowN <- results(dds_lowN, contrast = c(0,-1,0,1), alpha=0.05)
summary(res_W10vsW02_lowN)
degs_W10vsW02_lowN <- subset(res_W10vsW02_lowN, padj < 0.05)
dim(degs_W10vsW02_lowN)
degs_W10vsW02_lowN <- degs_W10vsW02_lowN[order(degs_W10vsW02_lowN$padj),]

# res_W05vsW10_lowN <- results(dds_lowN, contrast = c(0,0,1,-1), alpha=0.05)
# summary(res_W05vsW10_lowN)
# degs_W05vsW10_lowN <- subset(res_W05vsW10_lowN, padj < 0.05)
# dim(degs_W05vsW10_lowN)
# degs_W05vsW10_lowN <- degs_W05vsW10_lowN[order(degs_W05vsW10_lowN$padj),]
# 
# res_W01vsW02_lowN <- results(dds_lowN, contrast = c(1,-1,0,0), alpha=0.05)
# summary(res_W01vsW02_lowN)
# degs_W01vsW02_lowN <- subset(res_W01vsW02_lowN, padj < 0.05)
# dim(degs_W01vsW02_lowN)
# degs_W01vsW02_lowN <- degs_W01vsW02_lowN[order(degs_W01vsW02_lowN$padj),]

#Apply same filtering as above
filtered_degs_lowN <- rownames(degs_lowN) %>%
  intersect(rownames(degs_W05vsW01_lowN)) %>%
  intersect(rownames(degs_W05vsW02_lowN)) %>% 
  intersect(rownames(degs_W10vsW01_lowN)) %>%
  intersect(rownames(degs_W10vsW02_lowN)) #%>% 
  #setdiff(rownames(degs_W01vsW02_lowN)) %>% 
  #setdiff(rownames(degs_W05vsW10_lowN))

#filtered_degs_lowN_up <- rownames(degs_lowN[filtered_degs_lowN,][degs_lowN[filtered_degs_lowN,]$log2FoldChange>0,])
#filtered_degs_lowN_down <- rownames(degs_lowN[filtered_degs_lowN,][degs_lowN[filtered_degs_lowN,]$log2FoldChange<0,])

###Filter based on differential expression of homeologs
#1 Genes with all the homeologs being DE
all_homeologs_degs_lowN <- filtered_degs_lowN[sapply(filtered_degs_lowN,
                                                         function(deg) all_homeologs_appear(deg, homeologs,filtered_degs_lowN,F)) %>%
                                                    unlist()]

#1.1 Direction of fold-change (FC of same sign)
all_homeologs_degs_sameFC_lowN<- sapply(all_homeologs_degs_lowN, function(deg) all_homeologs_same_FC_sign(deg, homeologs, res_lowN))

#Save genes whose homeologs are DE in the same direction
#write.xlsx(homeolog_results_as_data_frame(all_homeologs_degs_lowN[all_homeologs_degs_sameFC_lowN], homeologs, res_lowN),
#           "filtered_results.xlsx",
#           sheetName = "lowN_all_homeologs_deg_sameFC",
#           row.names = F,
#           append = T)

#Save genes whose homeologs are DE in different directions
#write.xlsx(homeolog_results_as_data_frame(all_homeologs_degs_lowN[!all_homeologs_degs_sameFC_lowN], homeologs, res_lowN),
#           "filtered_results.xlsx",
#           sheetName = "lowN_all_homeologs_deg_diffFC",
#           row.names = F,
#           append = T)

#DEGs with no homeologs
singleton_degs_lowN <- filtered_degs_lowN[sapply(filtered_degs_lowN,
                                                   function(deg) is_singleton(deg))]

#write.xlsx(degs_lowN[singleton_degs_lowN, c("baseMean", "log2FoldChange", "padj")],
#           "filtered_results.xlsx",
#           sheetName = "lowN_singleton_degs",
#           row.names = T,
#           append = T)

#DEGs with more than one homeologs (but NOT all) DE
#Genes will be kept if:
# 1. None of the non-DE homeologs are expressed (were filtered out during processing)
# 2. The baseMean difference between the DEG with the lowest expression and the non-DE homeologs is more than one order of magnitude (>10x)
any_homeologs_degs_lowN <- filtered_degs_lowN[sapply(filtered_degs_lowN,
                                                       function(deg) any_homeologs_appear(deg, homeologs,filtered_degs_lowN,F)) %>%
                                                  unlist()]
multiple_homeologs_degs_lowN <- setdiff(any_homeologs_degs_lowN, all_homeologs_degs_lowN)
multiple_homeologs_degs_nondegs_low_lowN <- sapply(multiple_homeologs_degs_lowN,
                                                    function(gene){
                                                      gene_homeologs <- find_homeologs(gene, homeologs)
                                                      deg_homeologs <- gene_homeologs[gene_homeologs %in% multiple_homeologs_degs_lowN]
                                                      homeolog_results <- retrieve_homeolog_results(gene, homeologs, res_lowN)
                                                      min_deg_baseMean <- homeolog_results[deg_homeologs, "baseMean"] %>% min()
                                                      nondeg_results <- homeolog_results[! rownames(homeolog_results) %in% deg_homeologs,]
                                                      FC2DEG <- min_deg_baseMean/nondeg_results[, "baseMean"]
                                                      #print(gene)
                                                      if(FC2DEG %>% is.na() %>% sum()==length(gene_homeologs)-length(deg_homeologs) || all(FC2DEG > 10, na.rm = T) ) return(T) else return(F)})
multiple_homeologs_degs_nondegs_low_lowN <- multiple_homeologs_degs_lowN[multiple_homeologs_degs_nondegs_low_lowN]

write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs_nondegs_low_lowN, homeologs, res_lowN),
           "filtered_results.xlsx",
           sheetName = "lowN_multiple_homeologs_deg",
           row.names = F,
           append = T)



#DEGs with no DE homeologs but whose expression is much lower than that of the homeologs
#Genes will be kept if:
# 1. None of the other homeologs are expressed (were filtered out during processing)
# 2. The baseMean difference between the DEG and the homeologs is more than one order of magnitude (>10x)
no_homeologs_deg_lowN <- sapply(filtered_degs_lowN,
                                 function(deg) any_homeologs_appear(deg, homeologs,filtered_degs_lowN,T)) %>% unlist()
no_homeologs_deg_lowN <- filtered_degs_lowN[!no_homeologs_deg_lowN]

one_deg_low_homeologs_lowN <- sapply(no_homeologs_deg_lowN,
                                      function(gene){
                                        homeolog_results <- retrieve_homeolog_results(gene, homeologs, res_lowN) 
                                        FC2DEG <- homeolog_results[gene, "baseMean"]/homeolog_results[, "baseMean"]
                                        if(FC2DEG %>% is.na() %>% sum()==length(FC2DEG)-1 || all(FC2DEG > 10, na.rm = T) ) return(T) else return(F)})
one_deg_low_homeologs_lowN <- no_homeologs_deg_lowN[one_deg_low_homeologs_lowN]

write.xlsx(homeolog_results_as_data_frame(one_deg_low_homeologs_lowN, homeologs, res_lowN),
           "filtered_results.xlsx",
           sheetName = "lowN_one_deg_low_homeologs",
           row.names = F,
           append = T)


####################################################################
###########################FILTERING################################
####################################################################
source(file.path("//salt","wheat_rnaseq", "NRgeneV1", "ensembl_homeologs_utils.R"))

#Change as needed depending on what degs will be filtered
filtered_degs <- filtered_degs_lowN
res <- res_lowN

#1 Extract DEGs with no homeologs, aka singletons
singleton_degs <- filtered_degs[sapply(filtered_degs,
                                       function(deg) is_singleton(deg))]

#Remove singletons from the DEG list for further analysis
filtered_degs <- filtered_degs[! filtered_degs %in% singleton_degs]

#2 Genes with all the homeologs being DE
all_homeologs_degs <- filtered_degs[sapply(filtered_degs,
                                           function(deg){
                                             #gene_homeologs <- find_homeologs(deg)
                                             all_homeologs_appear(deg, filtered_degs,F)
                                           }) %>%
                                          unlist()]

#2.1 Direction of fold-change (FC of same sign)
all_homeologs_degs_sameFC<- sapply(all_homeologs_degs,
                                   function(deg){
                                     gene_homeologs <- find_homeologs(deg)
                                     genes_same_FC_sign(gene_homeologs,  res)
                                    } )

#Remove DEGS with all homeologs DE from the list for further analysis
filtered_degs <- filtered_degs[! filtered_degs %in% all_homeologs_degs]

#DEGs with more than one homeologs (but NOT all) DE
#DEGs are kept if:
# 1. None of the non-DE homeologs are expressed (were filtered out during processing)
# 2. The baseMean difference between the DEG with the lowest expression and the non-DE homeologs is more than one order of magnitude (>10x)

multiple_homeologs_degs <- filtered_degs[sapply(filtered_degs,
                                               function(deg) any_homeologs_appear(deg, filtered_degs,F)) %>%
                                          unlist()]

#Check the fold-change of the DE homeologs
multiple_homeologs_degs_sameFC<- sapply(multiple_homeologs_degs,
                                        function(deg) {
                                            gene_homeologs <- find_homeologs(deg)
                                            deg_homeologs <- gene_homeologs[gene_homeologs %in% multiple_homeologs_degs]
                                            genes_same_FC_sign(deg_homeologs,  res)})

multiple_homeologs_degs_nondegs_low <- sapply(multiple_homeologs_degs,
                                                function(gene){
                                                  gene_homeologs <- find_homeologs(gene)
                                                  deg_homeologs <- gene_homeologs[gene_homeologs %in% multiple_homeologs_degs]
                                                  #deg_homeologs <- c(gene, gene_homeologs[gene_homeologs %in% any_homeologs_degs_N])
                                                  homeolog_results <- retrieve_gene_results(gene_homeologs, res)
                                                  min_deg_baseMean <- homeolog_results[deg_homeologs, "baseMean"] %>% min()
                                                  nondeg_results <- homeolog_results[! rownames(homeolog_results) %in% deg_homeologs,]
                                                  FC2DEG <- min_deg_baseMean/nondeg_results[, "baseMean"]
                                                  #print(gene)
                                                  if(FC2DEG %>% is.na() %>% sum()==length(gene_homeologs)-length(deg_homeologs) || all(FC2DEG > 10, na.rm = T) ) return(T) else return(F)})

multiple_homeologs_degs_nondegs_low <- multiple_homeologs_degs[multiple_homeologs_degs_nondegs_low]

#Remove DEGS with multiple homeologs DE from the list for further analysis
#The degs left are those with no homeologs DE
filtered_degs <- filtered_degs[! filtered_degs %in% multiple_homeologs_degs]

#DEGs with no DE homeologs but whose expression is much higher than that of the homeologs
#Genes will be kept if:
# 1. None of the other homeologs are expressed (were filtered out during processing)
# 2. The baseMean difference between the DEG and the homeologs is more than one order of magnitude (>10x)

one_deg_low_homeologs <- sapply(filtered_degs,
                                  function(gene){
                                    gene_homeologs <- find_homeologs(gene)
                                    homeolog_results <- retrieve_gene_results(gene_homeologs, res) 
                                    FC2DEG <- homeolog_results[gene, "baseMean"]/homeolog_results[, "baseMean"]
                                    if(FC2DEG %>% is.na() %>% sum()==length(FC2DEG)-1 || all(FC2DEG > 10, na.rm = T) ) return(T) else return(F)})

one_deg_low_homeologs <- filtered_degs[one_deg_low_homeologs]


#Save all the results in an excel file
options(java.parameters = "- Xmx1024m")
library(xlsx)
#file_name <- "HIvLO_N_filtered_degs.xlsx" #Change as required
file_name <- "lowN_W01W02vW05W10_filtered_degs.xlsx" #Change as required

# 1. Singletons
if(length(singleton_degs)>0){
  write.xlsx(res[singleton_degs, c("baseMean", "log2FoldChange", "padj")],
             file_name,
             sheetName = "singleton_degs",
             row.names = T,
             append = F)
}


# 2. All homeologs DEG
# These genes are separated into two sheets depending on whether all the homeologs have the same FC sign

# 2.1 All homeologs DE and same FC sign
if(length(all_homeologs_degs[all_homeologs_degs_sameFC]) > 0){
  write.xlsx(homeolog_results_as_data_frame(all_homeologs_degs[all_homeologs_degs_sameFC], res),
             file_name,
             sheetName = "all_homeologs_deg_sameFC",
             row.names = F,
             append = T)
}


# 2.2 All homeologs DE but different FC sign
if(length(all_homeologs_degs[!all_homeologs_degs_sameFC]) > 0){
  write.xlsx(homeolog_results_as_data_frame(all_homeologs_degs[!all_homeologs_degs_sameFC], res),
           file_name,
           sheetName = "all_homeologs_deg_diffFC",
           row.names = F,
           append = T)
}

# 3. Multiple homeologs DE
# These genes are divided into three groups. First, all the genes with multiple homeologs DE are
# separated into two sheets depending on whether all the DE homeologs have the same FC sign. Then,
# the list genes with multiple homeologs is subseted to extract the cases in which the DE homeologs
# are more highly expressed than the non-DE homeologs

# 3.1 Genes with multiple homeologs DE and same FC
if(length(multiple_homeologs_degs[multiple_homeologs_degs_sameFC]) > 0){
  write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs[multiple_homeologs_degs_sameFC], res),
           file_name,
           sheetName = "multiple_homeologs_deg_sameFC",
           row.names = F,
           append = T)
}

# 3.2 Genes with multiple homeologs DE but different FC
if(length(multiple_homeologs_degs[!multiple_homeologs_degs_sameFC]) > 0){
  write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs[!multiple_homeologs_degs_sameFC], res),
             file_name,
             sheetName = "multiple_homeologs_deg_diffFC",
             row.names = F,
             append = T)
}

# 3.3 Genes with multiple homeologs DE but different FC
if(length(multiple_homeologs_degs_nondegs_low) > 0){
  write.xlsx(homeolog_results_as_data_frame(multiple_homeologs_degs_nondegs_low, res),
             file_name,
             sheetName = "multiple_homeologs_deg_lowNonDEGs",
             row.names = F,
             append = T)
}

# 4 Genes with a single homeologs DE
# First all the genes are saved and the a subdivision of genes where the DEG homeologs is more
# highly expressed than the non-DE homeologs

# 4.1 All the degs with no homeologs DE
if(length(filtered_degs)>0){
  write.xlsx(homeolog_results_as_data_frame(filtered_degs, res),
             file_name,
             sheetName = "one_homeolog_deg",
             row.names = F,
             append = T)
}

# 4.2 DEGs with no homeologs DE where the DEG homeolog is more highly expressed than the non-DEG
if(length(one_deg_low_homeologs)>0){
  write.xlsx(homeolog_results_as_data_frame(one_deg_low_homeologs, res),
             file_name,
             sheetName = "one_homeolog_deg_lowNonDEGs",
             row.names = F,
             append = T)
}


####################################################################
####################################################################
#W01 vs Paragon DEGs
res_W01 <- results(dds, name="variety_W01_vs_Paragon", alpha=0.05)
summary(res_W01)
degs_W01 <- subset(res_W01, padj < 0.05)
dim(degs_W01)
degs_W01 <- degs_W01[order(degs_W01$padj),]

#W02 vs Paragon DEGs
res_W02 <- results(dds, name="variety_W02_vs_Paragon", alpha=0.05)
summary(res_W02)
degs_W02 <- subset(res_W02, padj < 0.05)
dim(degs_W02)
degs_W02 <- degs_W02[order(degs_W02$padj),]

#W05 vs Paragon DEGs
res_W05 <- results(dds, name="variety_W05_vs_Paragon", alpha=0.05)
summary(res_W05)
degs_W05 <- subset(res_W05, padj < 0.05)
dim(degs_W05)
degs_W05 <- degs_W05[order(degs_W05$padj),]

#W10 vs Paragon DEGs
res_W10 <- results(dds, name="variety_W10_vs_Paragon", alpha=0.05)
summary(res_W10)
degs_W10 <- subset(res_W10, padj < 0.05)
dim(degs_W10)
degs_W10 <- degs_W10[order(degs_W10$padj),]

#Genes with upward trend
upward <- intersect(rownames(degs_W01[degs_W01$log2FoldChange<0,]), 
                    rownames(degs_W02[degs_W02$log2FoldChange<0,])) %>%
  intersect(rownames(degs_W05[degs_W05$log2FoldChange>0,])) %>%
  intersect(rownames(degs_W10[degs_W10$log2FoldChange>0,])) %>%
  intersect(rownames(degs_N))
length(upward)

#Genes with downward trend
downward <- intersect(rownames(degs_W01[degs_W01$log2FoldChange>0,]), 
                      rownames(degs_W02[degs_W02$log2FoldChange>0,])) %>%
  intersect(rownames(degs_W05[degs_W05$log2FoldChange<0,])) %>%
  intersect(rownames(degs_W10[degs_W10$log2FoldChange<0,])) %>%
  intersect(rownames(degs_N))
length(downward)

library(ggplot2)
#gene_id <- "TraesCS7B02G204400"
#sample_factors <- data.frame(variety=factor(sampleInfo$variety, levels(sampleInfo$variety)[c(2:3,1,4:5)]),
#                              condition = sampleInfo$condition,
#                              row.names = rownames(sampleInfo))
#dds <- dds_highN
#samples_factors <- highNInfo

##Functions for plotting the expression of the homeologs in the same graph
# When facet=F (default) all the homeologs will be plotted in the same graph. This is thought to be used when only
#one N level is plotted.
# When facet=T the homeologs will go in adjacent plots with high N and low N in the same graph.
plotHomoeologsExpression <- function(gene_vector, dds, sampleInfo, facet=F){
  
  #Assign a dummy expression value of 1 to those genes that do not appear
  #in dds, e.g. filtered due to low expression
  gene_in_results <- gene_vector %in% rownames(dds)
  NA_results <- data.frame(matrix(1, ncol=length(gene_vector[!gene_in_results]), nrow = nrow(sampleInfo)))
  colnames(NA_results) <- gene_vector[!gene_in_results]
  
  homeolog_expression <- cbind(sampleInfo, sapply(gene_vector[gene_in_results],
                                                  function(gene) counts(dds, normalized=T)[gene,]+0.5),
                              NA_results)
  if(facet){
    homeolog_expression <- homeolog_expression %>% gather(homeolog, norm_counts, -c(variety, condition, group))
    
    #Reorder varieties so that Paragon goes in the middle
    homeolog_expression$variety <- factor(homeolog_expression$variety, levels(homeolog_expression$variety)[c(2:3,1,4:5)])
    
    #Change the grouping variable to a combination of homeolog and condition
    homeolog_expression$group <- factor(paste(homeolog_expression$homeolog, homeolog_expression$condition, sep = "_"))
    
    x11()
    ggplot(homeolog_expression, aes(x=variety, y=norm_counts, group=group, col=condition)) +
      geom_point() +
      stat_summary(fun.y = "mean", geom = "line", size=1) +
      expand_limits(y=1) +
      scale_y_log10() +
      scale_colour_manual(values = c("blue", "red")) +
      facet_grid(~homeolog)
  }else{
    homeolog_expression <- homeolog_expression %>% gather(homeolog, norm_counts, -variety)
    x11()
    ggplot(homeolog_expression, aes(x=variety, y=norm_counts, group=homeolog, col=homeolog)) +
      geom_point() +
      stat_summary(fun.y = "mean", geom = "line", size=1) +
      expand_limits(y=1) +
      scale_y_log10()
  }
}

##Old version of functions

#Not flexible enough/inefficient - Deprecated
# plot2HomoeologsExpression <- function(title, geneA, geneB){
#   expression_plot_A <- plotExpression(geneA)
#   expression_plot_B <- plotExpression(geneB)
#   library("gridExtra")
#   x11()
#   grid.arrange(expression_plot_A, expression_plot_B, nrow=1, ncol=2, top = title)
# }

#Not flexible enough/inefficient - Deprecated
# plot3HomoeologsExpression <- function(title, geneA, geneB, geneD, dds, sample_factors){
#   sample_factors$norm_counts <- counts(dds, normalized=T)[geneA,]
#   expression_plot_A <- plotExpression(geneA, sample_factors)
#   sample_factors$norm_counts <- counts(dds, normalized=T)[geneB,]
#   expression_plot_B <- plotExpression(geneB, sample_factors)
#   sample_factors$norm_counts <- counts(dds, normalized=T)[geneD,]
#   expression_plot_D <- plotExpression(geneD, sample_factors)
#   library("gridExtra")
#   x11()
#   grid.arrange(expression_plot_A, expression_plot_B, expression_plot_D, nrow=1, ncol=3, top = title)
# }

#Too complex/inefficient - Deprecated
# deprecated_plotHomoeologsExpression <- function(gene_vector, dds, sample_factors){
#   if(length(gene_vector)<=1){
#     sample_factors$norm_counts <- counts(dds, normalized=T)[gene_vector,]
#     x11()
#     plotExpression(gene_vector, sample_factors)
#   }else{
#     expression_plots <- lapply(gene_vector,
#                                function(gene){ sample_factors$norm_counts <- counts(dds, normalized=T)[gene,]
#                                deprecated_plotExpression(gene, sample_factors)})
# 
#     #Calculate the number of rows and columns
#     ##If there are less than 4 genes plot them in one row
#     ##Otherwise put 3 plots per row
#     if(length(gene_vector)<=4){
#       rows=1
#       cols=length(gene_vector)
#     }else{
#       rows=ceiling(length(gene_vector)/3)
#       cols=3
#     }
#     x11();do.call(grid.arrange, c(expression_plots, nrow=rows, ncol=cols))
#   }
# }
# 
# #Too complex/inefficient - Deprecated
# deprecated_plotExpression <- function(gene_id, sample_factors){
#   #sample_factors$norm_counts <- counts(dds, normalized=T)[gene_id,]
#   #expression_plot <- ggplot(sample_factors, aes(x=variety, colour=condition)) +
#   expression_plot <- ggplot(sample_factors, aes(x=variety)) +
#     geom_point(aes(y=norm_counts)) +
#     scale_color_manual(values = c("blue", "red")) +
#     ggtitle(gene_id) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     expand_limits(y=0)
#   return(expression_plot)
# }