library(DESeq2)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(tidyverse)
library(GGally)

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load count matrixes
Bwa <- read.table("./counts/counts_bwa.txt", header=T, skip=1, row.names=1)[,6:11]
Hisat2Old <- read.table("./counts/counts_hisat2_old.txt", header=T, skip=1, row.names=1)[,6:11]
Hisat2 <- read.table("./counts/counts_hisat2.txt", header=T, skip=1, row.names=1)[,6:11]
Star   <- read.table("./counts/counts_star.txt", header=T, skip=1, row.names=1)[,6:11]
Tophat2 <- read.table("./counts/counts_tophat2.txt", header=T, skip=1, row.names=1)[,6:11]
Samples <- c("p15_s1","p15_s2","p15_s3","p5_s1","p5_s2","p5_s3")

colnames(Bwa) <- Samples
colnames(Hisat2Old)  <- Samples
colnames(Hisat2)  <- Samples
colnames(Star)  <- Samples
colnames(Tophat2) <- Samples

# groups
ColData <- data.frame(
  rownames = Samples,
  group = factor(c(rep('p5', 3), rep('p15', 3)))
)

# Functions
run_DESeq2 <- function(counts, coldata) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ group
  )
  dds <- dds[rowSums(counts(dds)) > 6, ]
  dds <- DESeq(dds)
  res <- results(dds, contrast = c('group', 'p15', 'p5'))
  return(res)
}

get_DGE <- function(res) {
  res <- as.data.frame(res)
  res <- na.omit(res)
  dge <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
  return(dge)
}

# DESeq2
ResBwa <- run_DESeq2(Bwa, ColData)
ResHisat2Old <- run_DESeq2(Hisat2Old, ColData)
ResHisat2 <- run_DESeq2(Hisat2, ColData)
ResStar <- run_DESeq2(Star, ColData)
ResTophat2 <- run_DESeq2(Tophat2, ColData)

# DGE
DgeBwa <- get_DGE(ResBwa)
DgeHisat2Old <- get_DGE(ResHisat2Old)
DgeHisat2 <- get_DGE(ResHisat2)
DgeStar <- get_DGE(ResStar)
DgeTophat2 <- get_DGE(ResTophat2)

# Venn diagram
venn.diagram(
  x = list(BWA = DgeBwa |> filter(log2FoldChange > 0) |> rownames(),
           HISAT2 = DgeHisat2 |> filter(log2FoldChange > 0) |> rownames(),
           HiSAT2NOSS = DgeHisat2Old |> filter(log2FoldChange > 0) |> rownames(),
           STAR = DgeStar |> filter(log2FoldChange > 0) |> rownames(),
           TopHat2 = DgeTophat2 |> filter(log2FoldChange > 0) |> rownames()
           ),
  filename = './VenDia.png',
  fill = c("tomato","steelblue","lightgreen","gold", 'purple3'), width = 3600, height = 3600
)


# DGE Spearman coef.
ResAll <-list(BWA = ResBwa, HISAT2NOSS = ResHisat2Old, HISAT2 = ResHisat2, STAR = ResStar, TopHat2 = ResTophat2) |>
  map(~as.data.frame(.x) |>
        rownames_to_column('gene') |>
        dplyr::select(gene, log2FoldChange) |>
        drop_na()) |>
  reduce(full_join, by = 'gene') |>
  setNames(c('gene', 'BWA', 'HISAT2NOSS', 'HISAT2', 'STAR', 'TopHat2')) |>
  drop_na() |>
  column_to_rownames('gene')

my_upper <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  r <- cor(x, y, method = "spearman")
  ggally_text(
    label = round(r, 3),
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    size = 1 + r * r * r * r * 16
  ) + theme(panel.background = element_rect(fill = "white"))
}

ggpairs(ResAll, upper = list(continuous = my_upper),
        lower = list(continuous = wrap("points", alpha = 0.05, size = 0.2)),
        diag  = list(continuous = wrap("densityDiag"))
) +
  theme_bw() +
  labs(title = "Pairwise log2FC Comparison (Spearman)")


##########
### GO ###
##########
ENSEMBL2ENTREZ <- function(deg_df) {
  genes_clean <- gsub("\\..*", "", rownames(deg_df))
  bitr(genes_clean,
       fromType = "ENSEMBL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db) |>
    pull(ENTREZID)
}

GOBwa <- enrichGO(gene = ENSEMBL2ENTREZ(DgeBwa |> filter(log2FoldChange > 0)), OrgDb = org.Hs.eg.db,
                  ont = 'BP', pAdjustMethod = 'BH', readable = T)

GOHisat2Old <- enrichGO(gene = ENSEMBL2ENTREZ(DgeHisat2Old |> filter(log2FoldChange > 0)), OrgDb = org.Hs.eg.db,
                  ont = 'BP', pAdjustMethod = 'BH', readable = T)

GOHisat2 <- enrichGO(gene = ENSEMBL2ENTREZ(DgeHisat2 |> filter(log2FoldChange > 0)), OrgDb = org.Hs.eg.db,
                  ont = 'BP', pAdjustMethod = 'BH', readable = T)

GOStar <- enrichGO(gene = ENSEMBL2ENTREZ(DgeStar |> filter(log2FoldChange > 0)), OrgDb = org.Hs.eg.db,
                  ont = 'BP', pAdjustMethod = 'BH', readable = T)

GOTophat2 <- enrichGO(gene = ENSEMBL2ENTREZ(DgeTophat2 |> filter(log2FoldChange > 0)), OrgDb = org.Hs.eg.db,
                  ont = 'BP', pAdjustMethod = 'BH', readable = T)

dotplot(GOBwa,     showCategory = 20, title = "GO BP - BWA")
dotplot(GOHisat2Old,  showCategory = 20, title = "GO BP - HISAT2NOSS")
dotplot(GOHisat2,    showCategory = 20, title = "GO BP - HISAT2")
dotplot(GOStar, showCategory = 20, title = "GO BP - STAR")
dotplot(GOTophat2, showCategory = 20, title = "GO BP - TopHat2")

# Union
deg_list <- list(
  BWA        = ENSEMBL2ENTREZ(DgeBwa |> filter(log2FoldChange > 0)),
  HISAT2NOSS = ENSEMBL2ENTREZ(DgeHisat2Old |> filter(log2FoldChange > 0)),
  HISAT2     = ENSEMBL2ENTREZ(DgeHisat2 |> filter(log2FoldChange > 0)),
  STAR       = ENSEMBL2ENTREZ(DgeStar |> filter(log2FoldChange > 0)),
  TopHat2    = ENSEMBL2ENTREZ(DgeTophat2 |> filter(log2FoldChange > 0))
)

ck <- compareCluster(deg_list, fun = "enrichGO",
                     OrgDb = org.Hs.eg.db,
                     ont = "BP", pAdjustMethod = "BH", readable = TRUE)

dotplot(ck, showCategory = 10) +
  theme_bw() +
  labs(title = "GO BP Comparison across Aligners")

venn.diagram(
  x = list(BWA = GOBwa@result |> filter(p.adjust < 0.05) |> _$ID,
           HISAT2 = GOHisat2@result |> filter(p.adjust < 0.05) |> _$ID,
           HISAT2NOSS = GOHisat2Old@result |> filter(p.adjust < 0.05) |> _$ID,
           STAR = GOStar@result |> filter(p.adjust < 0.05) |> _$ID,
           TopHat2 = GOTophat2@result |> filter(p.adjust < 0.05) |> _$ID
  ),
  filename = './VenDia_GOBP.png',
  fill = c("tomato","steelblue","lightgreen","gold", 'purple3'), width = 3600, height = 3600
)

save.image('./DGE_GO.RData')
