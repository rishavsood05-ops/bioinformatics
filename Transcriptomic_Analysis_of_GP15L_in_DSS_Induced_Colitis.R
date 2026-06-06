library(DESeq2)
library(tidyverse)

#Expression Data

counts <- read.csv("C:/Users/risha/Downloads/GSE327543_Fig3ab.csv.gz", row.names = 1, check.names = FALSE)

head(counts)
dim(counts)


#making Counts numeric

counts <- as.matrix(counts)

mode(counts) <- "numeric"

counts <- round(counts)
head(counts)

sample_info <- data.frame(
  sample = colnames(counts),
  condition = c(
    "Colon_WT_DSS",
    "Colon_WT_DSS",
    "Colon_WT_DSS",
    "Colon_GPR15L_KO_DSS",
    "Colon_GPR15L_KO_DSS",
    "Colon_GPR15L_KO_DSS",
    "LPMCs_WT_DSS",
    "LPMCs_WT_DSS",
    "LPMCs_WT_DSS",
    "LPMCs_GPR15L_KO_DSS",
    "LPMCs_GPR15L_KO_DSS",
    "LPMCs_GPR15L_KO_DSS"
  )
)

rownames(sample_info) <- sample_info$sample
sample_info


#error checking
colnames(counts)
rownames(sample_info)

head(counts)
sum(is.na(counts))

dim(counts)

head(counts[,1:4])

colnames(counts)

sample_info

str(counts)
#run DESeq2

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Colon_GPR15L_KO_DSS", "Colon_WT_DSS")
)

res_df <- as.data.frame(res)
head(res_df)

head(res_df, 20)

#sort by significance

res_df <- as.data.frame(res)

res_df <- res_df[order(res_df$padj), ]

head(res_df, 20)

#no. of significant genes

sum(res_df$padj < 0.05, na.rm = TRUE)

# no. of adjusted p. values

sum(!is.na(res_df$padj))

#top significant genes

sig_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05)

head(sig_genes, 20)

#Creating the volcano plot

library(ggplot2)

ggplot(res_df,
       aes(x = log2FoldChange,
           y = -log10(pvalue))) +
  geom_point(alpha = 0.5) +
  theme_minimal()

nrow(res_df)

summary(res_df$padj)

sum(res_df$padj < 0.05, na.rm = TRUE)

#Redoing the volcano plot

volcano_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    Significant = ifelse(
      padj < 0.05 & abs(log2FoldChange) > 1,
      "Yes",
      "No"
    )
  )

ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = Significant)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

#finding most significant genes
library(dplyr)

top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

head(top_genes, 20)

#creating the heatmap

library(dplyr)
library(pheatmap)

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

top30 <- res_df$gene[1:30]

vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)

heatmap_data <- norm_counts[top30, ]

heatmap_data <- heatmap_data[
  complete.cases(heatmap_data),
]

heatmap_data <- heatmap_data[
  apply(heatmap_data, 1, var) > 0,
]

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = TRUE,
  main = "Top 30 Differentially Expressed Genes"
)


#download org.Mm.eg.db and AnnotationDbi
BiocManager::install("org.Mm.eg.db")
BiocManager::install("AnnotationDbi")

#convert Ensembles to actual gene names

library(org.Mm.eg.db)
library(AnnotationDbi)

top_genes$symbol <- mapIds(
  org.Mm.eg.db,
  keys = top_genes$gene,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

head(top_genes[, c("gene","symbol","log2FoldChange","padj")], 20)



resultsNames(dds)

#Separating Colon samples from LPMC

res_colon <- results(
  dds,
  name = "condition_Colon_WT_DSS_vs_Colon_GPR15L_KO_DSS"
)

res_colon_df <- as.data.frame(res_colon) %>%
  rownames_to_column("gene") %>%
  arrange(padj)
res_colon_KO_vs_WT <- results(
  dds,
  contrast = c("condition", "Colon_GPR15L_KO_DSS", "Colon_WT_DSS")
)

res_colon_KO_vs_WT_df <- as.data.frame(res_colon_KO_vs_WT) %>%
  rownames_to_column("gene") %>%
  arrange(padj)
head(res_colon_KO_vs_WT_df)

summary(res_colon_KO_vs_WT)

head(
  res_colon_KO_vs_WT_df[
    order(res_colon_KO_vs_WT_df$padj),
  ],
  20
)
levels(dds$condition)

head(res_colon_KO_vs_WT_df)

#Volcano Plot for Colon Samples

library(ggplot2)
library(dplyr)
library(tibble)

res_colon <- results(
  dds,
  contrast = c("condition", "Colon_GPR15L_KO_DSS", "Colon_WT_DSS")
)

res_colon_df <- as.data.frame(res_colon) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj)) %>%
  mutate(
    Significant = ifelse(
      padj < 0.05 & abs(log2FoldChange) > 1,
      "Yes",
      "No"
    )
  )

volcano_colon <- ggplot(
  res_colon_df,
  aes(x = log2FoldChange, y = -log10(padj), color = Significant)
) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Colon GPR15L-KO DSS vs Colon WT DSS",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

print(volcano_colon)

# Colon Heatmap

library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)

top30_colon <- res_colon_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:30) %>%
  pull(gene)

colon_samples <- c(
  "Colon_WT_DSS_1",
  "Colon_WT_DSS_2",
  "Colon_WT_DSS_3",
  "Colon_GPR15L-KO_DSS_1",
  "Colon_GPR15L-KO_DSS_2",
  "Colon_GPR15L-KO_DSS_3"
)

heatmap_colon <- norm_counts[top30_colon, colon_samples]

heatmap_colon <- heatmap_colon[
  complete.cases(heatmap_colon),
]

heatmap_colon <- heatmap_colon[
  apply(heatmap_colon, 1, var) > 0,
]

annotation_colon <- data.frame(
  condition = c(
    "WT",
    "WT",
    "WT",
    "GPR15L_KO",
    "GPR15L_KO",
    "GPR15L_KO"
  )
)

rownames(annotation_colon) <- colon_samples

pheatmap(
  heatmap_colon,
  scale = "row",
  annotation_col = annotation_colon,
  show_rownames = TRUE,
  main = "Top 30 DE Genes: Colon KO vs WT"
)

#LPMC Volcano Plot

res_LPMC <- results(
  dds,
  contrast = c("condition", "LPMCs_GPR15L_KO_DSS", "LPMCs_WT_DSS")
)

res_LPMC_df <- as.data.frame(res_LPMC) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

volcano_LPMC <- res_LPMC_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    Significant = ifelse(
      padj < 0.05 & abs(log2FoldChange) > 1,
      "Yes",
      "No"
    )
  )

ggplot(volcano_LPMC,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = Significant)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "LPMCs GPR15L-KO DSS vs LPMCs WT DSS",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )

#LPMC Heatmap

top30_LPMC <- res_LPMC_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:30) %>%
  pull(gene)

LPMC_samples <- c(
  "LPMCs_WT_DSS_1",
  "LPMCs_WT_DSS_2",
  "LPMCs_WT_DSS_3",
  "LPMCs_GPR15L-KO_DSS_1",
  "LPMCs_GPR15L-KO_DSS_2",
  "LPMCs_GPR15L-KO_DSS_3"
)

heatmap_LPMC <- norm_counts[top30_LPMC, LPMC_samples]

heatmap_LPMC <- heatmap_LPMC[
  complete.cases(heatmap_LPMC),
]

heatmap_LPMC <- heatmap_LPMC[
  apply(heatmap_LPMC, 1, var) > 0,
]

annotation_LPMC <- data.frame(
  condition = c(
    "WT",
    "WT",
    "WT",
    "GPR15L_KO",
    "GPR15L_KO",
    "GPR15L_KO"
  )
)

rownames(annotation_LPMC) <- LPMC_samples

pheatmap(
  heatmap_LPMC,
  scale = "row",
  annotation_col = annotation_LPMC,
  show_rownames = TRUE,
  main = "Top 30 DE Genes: LPMC KO vs WT"
)


#Enrichment Pathway for Colon

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)


sig_colon <- res_colon_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) > 1
  )

nrow(sig_colon)



gene_conversion <- bitr(
  sig_colon$gene,
  fromType = "ENSEMBL",
  toType = c("ENTREZID","SYMBOL"),
  OrgDb = org.Mm.eg.db
)

head(gene_conversion)

go_colon <- enrichGO(
  gene = gene_conversion$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

head(go_colon)

as.data.frame(go_colon)


#Barplot for Colon Results

barplot(
  go_colon,
  showCategory = 15,
  title = "GO Biological Processes"
)

dotplot(
  go_colon,
  showCategory = 15
)


#KEGG enrichment pathway for Colon

kegg_colon <- enrichKEGG(
  gene = gene_conversion$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)

head(kegg_colon)


as.data.frame(kegg_colon)

dotplot(
  kegg_colon,
  showCategory = 15
)

#Enrichment Pathway for LPMC

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)


sig_LPMC <- res_LPMC_df %>%
  filter(
    !is.na(padj),
    padj < 0.05,
    abs(log2FoldChange) > 1
  )

gene_conversion_LPMC <- bitr(
  sig_LPMC$gene,
  fromType = "ENSEMBL",
  toType = c("ENTREZID","SYMBOL"),
  OrgDb = org.Mm.eg.db
)

go_LPMC <- enrichGO(
  gene = gene_conversion_LPMC$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP"
)

head(go_LPMC)
as.data.frame(go_LPMC)

#Barplot for LPMC Results

barplot(
  go_LPMC,
  showCategory = 15,
  title = "GO Biological Processes"
)

dotplot(
  go_LPMC,
  showCategory = 15
)

#LPMC additional data

nrow(sig_LPMC)

nrow(gene_conversion_LPMC)

head(as.data.frame(go_LPMC))

head(as.data.frame(kegg_LPMC))

#Identifying Significant genes for LPMC

sig_LPMC <- res_LPMC_df %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

nrow(sig_LPMC)

res_LPMC_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(10)

top_LPMC <- res_LPMC_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:10)

mapIds(
  org.Mm.eg.db,
  keys = top_LPMC$gene,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
