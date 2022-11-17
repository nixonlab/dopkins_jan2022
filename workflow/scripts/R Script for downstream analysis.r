title: "Nick HERV Micro Proper"
output: html_document
---

```{r setup, include=FALSE}
library("knitr")
knitr::opts_chunk$set(warning = FALSE)
library("BiocManager")
library("readxl")
library("gplots")
library("ggplot2")
library("DESeq2")
library("edgeR")
library("cowplot")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("biomaRt")
library("org.Hs.eg.db")
library("AnnotationDbi")
```

## Importing TE Data
```{r}
HERV_Counts <- read.csv(file = "HERV_Counts.csv", row.names = 1, header = TRUE)

colnames(HERV_Counts) <- c('FLA_Lymph 1', 'FLA_Mono 1', 'LPS_Lymph 1', 'LPS_Mono 1', '0h_Lymph 1', '0h_Mono 1', '2h_Lymph 1', '2h_Mono 1', 'FLA_Lymph 2', 'FLA_Mono 2', 'LPS_Lymph 2', 'LPS_Mono 2', '0h_Lymph 2', '0h_Mono 2', '2h_Lymph 2', '2h_Mono 2', 'FLA_Lymph 3', 'FLA_Mono 3', 'LPS_Lymph 3', 'LPS_Mono 3', '0h_Lymph 3', '0h_Mono 3', '2h_Lymph 3', '2h_Mono 3')

HERVs <- as.data.frame(HERV_Counts) 
HERV_Counts <- c("FLA_Lymph 1", "FLA_Mono 1", "LPS_Lymph 1", "LPS_Mono 1", "0h_Lymph 1", "0h_Mono 1", "2h_Lymph 1", "2h_Mono 1", "FLA_Lymph 2", "FLA_Mono 2", "LPS_Lymph 2", "LPS_Mono 2", "0h_Lymph 2", "0h_Mono 2", "2h_Lymph 2", "2h_Mono 2", "FLA_Lymph 3", "FLA_Mono 3", "LPS_Lymph 3", "LPS_Mono 3", "0h_Lymph 3", "0h_Mono 3", "2h_Lymph 3", "2h_Mono 3")


```

##DESEQ counts
```{r}
HERVsDESEQ <- HERVs[!(row.names(HERVs) %in% "__no_feature"),]
HERVsMonoDESEQ <- HERVsDESEQ[c(2,4,6,8,10,12,14,16,18,20,22,24)]
HERVsLymphoDESEQ <- HERVsDESEQ[c(1,3,5,7,9,11,13,15,17,19,21,23)]
cutoff.count <- 5

##Minimum number of samples meeting the minimum count threshold
cutoff.samp1 <- floor(ncol(HERVsMonoDESEQ) * 0.20)
HERVsMonoDESEQ <- HERVsMonoDESEQ[rowSums(HERVsMonoDESEQ > cutoff.count) > cutoff.samp1, ]

cutoff.samp2 <- floor(ncol(HERVsLymphoDESEQ) * 0.20)
HERVsLymphoDESEQ <- HERVsLymphoDESEQ[rowSums(HERVsLymphoDESEQ > cutoff.count) > cutoff.samp2, ]

##MonoDESEQ and Volcano Plots

ddsfMono <- read.csv(file = "DESeq2KeyMono.csv")
HERVSCPMMonoMatrixDESEQ <- as.matrix(HERVsMonoDESEQ)
ddsHERVsMonoDESEQ <- DESeqDataSetFromMatrix(countData =  HERVSCPMMonoMatrixDESEQ, colData = ddsfMono, design = ~ Donor + Condition + Donor: Condition)
ddsHERVsMonoDESEQ <- DESeq(ddsHERVsMonoDESEQ)
resultsNames(ddsHERVsMonoDESEQ)
resultsMonoFLADESEQ <- results(ddsHERVsMonoDESEQ, name = "Condition_FLA_vs_CTRL")
resultsMonoLPSDESEQ <- results(ddsHERVsMonoDESEQ, name = "Condition_LPS_vs_CTRL")

pdf("FLA_Mono_Volcano_Plot1.pdf")
EnhancedVolcano::EnhancedVolcano(resultsMonoFLADESEQ, lab = rownames(resultsMonoFLADESEQ), x = 'log2FoldChange', y = 'pvalue', title = 'Flagella Stimulated Monocytes', pCutoff = 0.05, FCcutoff = 1, labSize = 2.5,  labFace = 'bold', xlim = c(-20, 20), ylim = c(0, 5))

dev.off()

pdf("LPS_Mono_Volcano_Plot1.pdf")
EnhancedVolcano::EnhancedVolcano(resultsMonoLPSDESEQ, lab = rownames(resultsMonoLPSDESEQ), x = 'log2FoldChange', y = 'pvalue', title = 'LPS Stimulated Monocytes', pCutoff = 0.05, FCcutoff = 1, labSize = 2.5,  labFace = 'bold', xlim = c(-20, 20), ylim = c(0, 5))

dev.off()



##LymphoDESEQ and Volcano Plots
ddsfLympho <- read.csv(file = "DESeq2KeyLympho.csv")
HERVSCPMLymphoMatrixDESEQ <- as.matrix(HERVsLymphoDESEQ)
ddsHERVsLymphoDESEQ <- DESeqDataSetFromMatrix(countData =  HERVSCPMLymphoMatrixDESEQ, colData = ddsfMono, design = ~ Donor + Condition + Donor: Condition)
ddsHERVsLymphoDESEQ <- DESeq(ddsHERVsLymphoDESEQ)
resultsNames(ddsHERVsLymphoDESEQ)
resultsLymphoFLADESEQ <- results(ddsHERVsLymphoDESEQ, name = "Condition_FLA_vs_CTRL")
resultsLymphoLPSDESEQ <- results(ddsHERVsLymphoDESEQ, name = "Condition_LPS_vs_CTRL")

pdf("FLA_Lympho_Volcano_Plot1.pdf")
EnhancedVolcano::EnhancedVolcano(resultsLymphoFLADESEQ, lab = rownames(resultsLymphoFLADESEQ), x = 'log2FoldChange', y = 'pvalue', title = 'Flagella Stimulated Lymphocytes', pCutoff = 0.05, FCcutoff = 1, labSize = 2.5,  labFace = 'bold', xlim = c(-20, 20), ylim = c(0, 5))

dev.off()

pdf("LPS_Lympho_Volcano_Plot1.pdf")
EnhancedVolcano::EnhancedVolcano(resultsLymphoLPSDESEQ, lab = rownames(resultsLymphoLPSDESEQ), x = 'log2FoldChange', y = 'pvalue', title = 'LPS Stimulated Lymphocytes', pCutoff = 0.05, FCcutoff = 1, labSize = 2.5,  labFace = 'bold',, xlim = c(-20, 20), ylim = c(0, 5))

dev.off()


##GTF Files for ENSG

```{r}
load("herv_annot.Rdata") 
retro.annot <- read.table("genes.tsv", sep='\t', header=T)
retro.annot$family <- sapply(strsplit(retro.annot$locus, '_'), '[[', 1)
row.names(retro.annot) <- retro.annot$locus
retro.annot$family <- annot.herv$family[match(retro.annot$locus, annot.herv$locus)]
retro.annot$group <- annot.herv$group[match(retro.annot$locus, annot.herv$locus)]
retro.annot$letter <- annot.herv$letter[match(retro.annot$locus, annot.herv$locus)]

# Add L1 as a factor, and change the family, group, and letter for the L1s
retro.annot$family = factor(retro.annot$family, levels=c(levels(retro.annot$family), "L1"))
retro.annot$family[is.na(retro.annot$family)] = "L1"
retro.annot$group = factor(retro.annot$group, levels=c(levels(retro.annot$group), "L1"))
retro.annot$group[is.na(retro.annot$group)] = "L1"
retro.annot$letter = factor(retro.annot$letter, levels=c(levels(retro.annot$letter), "L1"))
retro.annot$letter[is.na(retro.annot$letter)] = "L1"

# Add LINE1 and HERVs as classes 
retro.annot$class <- ifelse(grepl('^L1', rownames(retro.annot)), 'LINE1', 'HERV')
retro.annot$class <- factor(retro.annot$class)


load("gene_annot.Rdata")
##Create gene ID to gene name mapping file
gene_table <- gtf_df[!duplicated(gtf_df[,c(1,7)]), ] %>% dplyr::select('gene_id', 'gene_name', 'gene_type')
gene_table <- rbind(gene_table, data.frame(gene_id=retro.annot$locus, gene_name=retro.annot$locus, gene_type=retro.annot$class))
rownames(gene_table) <- gene_table$gene_id

gene_table$display <- gene_table$gene_name
gene_table[duplicated(gene_table$gene_name), 'display'] <- paste(gene_table[duplicated(gene_table$gene_name), 'display'], gene_table[duplicated(gene_table$gene_name), 'gene_id'], sep='|')


```


##Global Gene DESEQ and Volcano Plots

```{r}

Gene_Counts <- read.csv(file = "Gene_counts.csv", row.names = 1, header = TRUE)

colnames(Gene_Counts) <- c('FLA_Lymph 1', 'FLA_Mono 1', 'LPS_Lymph 1', 'LPS_Mono 1', '0h_Lymph 1', '0h_Mono 1', '2h_Lymph 1', '2h_Mono 1', 'FLA_Lymph 2', 'FLA_Mono 2', 'LPS_Lymph 2', 'LPS_Mono 2', '0h_Lymph 2', '0h_Mono 2', '2h_Lymph 2', '2h_Mono 2', 'FLA_Lymph 3', 'FLA_Mono 3', 'LPS_Lymph 3', 'LPS_Mono 3', '0h_Lymph 3', '0h_Mono 3', '2h_Lymph 3', '2h_Mono 3')

Genes <- as.data.frame(Gene_Counts) 
Gene_Counts <- c("FLA_Lymph 1", "FLA_Mono 1", "LPS_Lymph 1", "LPS_Mono 1", "0h_Lymph 1", "0h_Mono 1", "2h_Lymph 1", "2h_Mono 1", "FLA_Lymph 2", "FLA_Mono 2", "LPS_Lymph 2", "LPS_Mono 2", "0h_Lymph 2", "0h_Mono 2", "2h_Lymph 2", "2h_Mono 2", "FLA_Lymph 3", "FLA_Mono 3", "LPS_Lymph 3", "LPS_Mono 3", "0h_Lymph 3", "0h_Mono 3", "2h_Lymph 3", "2h_Mono 3")
# No group or design set. Assuming all samples belong to one group.


GenesDESEQ <- Genes[!(row.names(Genes) %in% "N_ambiguous"),]
GenesDESEQ <- Genes[!(row.names(Genes) %in% "N_noFeature"),]

GenesMonoDESEQ <- GenesDESEQ[c(2,4,6,8,10,12,14,16,18,20,22,24)]
GenesLymphoDESEQ <- GenesDESEQ[c(1,3,5,7,9,11,13,15,17,19,21,23)]
cutoff.count <- 5

# Minimum number of samples meeting the minimum count threshold

GenesMonoDESEQ <- GenesMonoDESEQ[rowSums(GenesMonoDESEQ > cutoff.count) > cutoff.samp1, ]
GenesLymphoDESEQ <- GenesLymphoDESEQ[rowSums(GenesLymphoDESEQ > cutoff.count) > cutoff.samp2, ]


##MonoDESEQ and Volcano Plots


GeneSCPMMonoMatrixDESEQ <- as.matrix(GenesMonoDESEQ)
ddsGenesMonoDESEQ <- DESeqDataSetFromMatrix(countData =  GeneSCPMMonoMatrixDESEQ, colData = ddsfMono, design = ~ Donor + Condition + Donor: Condition)
ddsGenesMonoDESEQ <- DESeq(ddsGenesMonoDESEQ)
resultsNames(ddsGenesMonoDESEQ)
resultsMonoFLADESEQGenes <- results(ddsGenesMonoDESEQ, name = "Condition_FLA_vs_CTRL")
resultsMonoLPSDESEQGenes <- results(ddsGenesMonoDESEQ, name = "Condition_LPS_vs_CTRL")



pdf("FLA_Mono_Volcano_Plot_Genes.pdf")
EnhancedVolcano::EnhancedVolcano(resultsMonoFLADESEQGenes, lab = gene_table[rownames(resultsMonoFLADESEQGenes), 'display'], x = 'log2FoldChange', y = 'pvalue', title = 'Flagella Stimulated Monocytes', pCutoff = 0.05, FCcutoff = 1, labSize = 1.5, xlim = c(-20, 20), ylim = c(0, 5))

dev.off()

pdf("LPS_Mono_Volcano_Plot_Genes.pdf")
EnhancedVolcano::EnhancedVolcano(resultsMonoLPSDESEQGenes, lab = gene_table[rownames(resultsMonoLPSDESEQGenes), 'display'], x = 'log2FoldChange', y = 'pvalue', title = 'LPS Stimulated Monocytes', pCutoff = 0.05, FCcutoff = 1, labSize = 1.5, xlim = c(-20, 20), ylim = c(0, 5))

dev.off()

##LymphoDESEQ and Volcano Plots


GeneSCPMLymphoMatrixDESEQ <- as.matrix(GenesLymphoDESEQ)
ddsGenesLymphoDESEQ <- DESeqDataSetFromMatrix(countData =  GeneSCPMLymphoMatrixDESEQ, colData = ddsfLympho, design = ~ Donor + Condition + Donor: Condition)
ddsGenesLymphoDESEQ <- DESeq(ddsGenesLymphoDESEQ)
resultsNames(ddsGenesLymphoDESEQ)
resultsLymphoFLADESEQGenes <- results(ddsGenesLymphoDESEQ, name = "Condition_FLA_vs_CTRL")
resultsLymphoLPSDESEQGenes <- results(ddsGenesLymphoDESEQ, name = "Condition_LPS_vs_CTRL")

pdf("FLA_Lympho_Volcano_Plot_Genes.pdf")
EnhancedVolcano::EnhancedVolcano(resultsLymphoFLADESEQGenes, lab = gene_table[rownames(resultsLymphoFLADESEQGenes), 'display'], x = 'log2FoldChange', y = 'pvalue', title = 'Flagella Stimulated Lymphocytes', pCutoff = 0.05, FCcutoff = 1, labSize = 1.5, xlim = c(-20, 20), ylim = c(0, 5))

dev.off()

pdf("LPS_Lympho_Volcano_Plot_Genes.pdf")
EnhancedVolcano::EnhancedVolcano(resultsLymphoLPSDESEQGenes, lab = gene_table[rownames(resultsLymphoLPSDESEQGenes), 'display'], x = 'log2FoldChange', y = 'pvalue', title = 'LPS Stimulated Lymphoocytes', pCutoff = 0.05, FCcutoff = 1, labSize = 1.5, xlim = c(-20, 20), ylim = c(0, 5))

dev.off()



##Export significant genes from DESEQ
```{r}
##lymphoFLA
resultsLymphoFLADESEQGenesignificant <- subset(resultsLymphoFLADESEQGenes, pvalue < 0.05)
resultsLymphoFLADESEQGenesignificant <- subset(resultsLymphoFLADESEQGenesignificant, log2FoldChange >= 1 | log2FoldChange <= -1)
resultsLymphoFLADESEQGenesignificant <- as.data.frame(resultsLymphoFLADESEQGenesignificant)
resultsLymphoFLADESEQGenesignificant$gene_names <- gene_table[rownames(resultsLymphoFLADESEQGenesignificant),"gene_name"]
write.csv(resultsLymphoFLADESEQGenesignificant,"C:/Users/17729/Documents/Nick HERV Micro/LymphoFLAGenesSig", row.names = TRUE)

##monoFLA

resultsMonoFLADESEQGenesignificant <- subset(resultsMonoFLADESEQGenes, pvalue < 0.05)
resultsMonoFLADESEQGenesignificant <- subset(resultsMonoFLADESEQGenesignificant, log2FoldChange >= 1 | log2FoldChange <= -1)
resultsMonoFLADESEQGenesignificant <- as.data.frame(resultsMonoFLADESEQGenesignificant)
resultsMonoFLADESEQGenesignificant$gene_names <- gene_table[rownames(resultsMonoFLADESEQGenesignificant),"gene_name"]
write.csv(resultsMonoFLADESEQGenesignificant,"C:/Users/17729/Documents/Nick HERV Micro/MonoFLAGenesSig", row.names = TRUE)

##lymphoLPS
resultsLymphoLPSDESEQGenesignificant <- subset(resultsLymphoLPSDESEQGenes, pvalue < 0.05)
resultsLymphoLPSDESEQGenesignificant <- subset(resultsLymphoLPSDESEQGenesignificant, log2FoldChange >= 1 | log2FoldChange <= -1)
resultsLymphoLPSDESEQGenesignificant <- as.data.frame(resultsLymphoLPSDESEQGenesignificant)
resultsLymphoLPSDESEQGenesignificant$gene_names <- gene_table[rownames(resultsLymphoLPSDESEQGenesignificant),"gene_name"]
write.csv(resultsLymphoLPSDESEQGenesignificant,"C:/Users/17729/Documents/Nick HERV Micro/LymphoLPSGenesSig", row.names = TRUE)

##monoLPS

resultsMonoLPSDESEQGenesignificant <- subset(resultsMonoLPSDESEQGenes, pvalue < 0.05)
resultsMonoLPSDESEQGenesignificant <- subset(resultsMonoLPSDESEQGenesignificant, log2FoldChange >= 1 | log2FoldChange <= -1)
resultsMonoLPSDESEQGenesignificant <- as.data.frame(resultsMonoLPSDESEQGenesignificant)
resultsMonoLPSDESEQGenesignificant$gene_names <- gene_table[rownames(resultsMonoLPSDESEQGenesignificant),"gene_name"]
write.csv(resultsMonoLPSDESEQGenesignificant,"C:/Users/17729/Documents/Nick HERV Micro/MonoLPSGenesSig", row.names = TRUE)



```
## Export CSV of significantly modified TEs

```{r}

resultslymphoFLAsignificant <- subset(resultsLymphoFLADESEQ, pvalue < 0.05)
resultslymphoFLAsignificant <- subset(resultslymphoFLAsignificant, log2FoldChange >= 1 | log2FoldChange <= -1)
resultslymphoFLAsignificant<- as.data.frame(resultslymphoFLAsignificant)
write.csv(resultslymphoFLAsignificant,"C:/Users/17729/Documents/Nick HERV Micro/TEsLymphFLA", row.names = TRUE)


resultslymphoLPSsignificant <- subset(resultsLymphoLPSDESEQ, pvalue < 0.05)
resultslymphoLPSsignificant <- subset(resultslymphoLPSsignificant, log2FoldChange >= 1 | log2FoldChange <= -1)
resultslymphoLPSsignificant <- as.data.frame(resultslymphoLPSsignificant)
write.csv(resultslymphoLPSsignificant,"C:/Users/17729/Documents/Nick HERV Micro/TEsLymphLPS", row.names = TRUE)

resultsmonoFLAsignificant <- subset(resultsMonoFLADESEQ, pvalue < 0.05)
resultsmonoFLAsignificant <- subset(resultsmonoFLAsignificant, log2FoldChange >= 1 | log2FoldChange <= -1)
resultsmonoFLAsignificant <- as.data.frame(resultsmonoFLAsignificant)
write.csv(resultsmonoFLAsignificant,"C:/Users/17729/Documents/Nick HERV Micro/TEsMonoFLA", row.names = TRUE)

resultsmonoLPSsignificant <- subset(resultsMonoLPSDESEQ, pvalue < 0.05)
resultsmonoLPSsignificant <- subset(resultsmonoLPSsignificant, log2FoldChange >= 1 | log2FoldChange <= -1)
resultsmonoLPSsignificant <- as.data.frame(resultsmonoLPSsignificant)
write.csv(resultsmonoLPSsignificant,"C:/Users/17729/Documents/Nick HERV Micro/TEsMonoLPS", row.names = TRUE)


```

