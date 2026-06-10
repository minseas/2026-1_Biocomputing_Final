library(DESeq2)
library(ggplot2)
library(tidyverse)

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


# combine
counts_list <- list(
  BWA = Bwa,
  HISAT2NOSS = Hisat2Old,
  HISAT2 = Hisat2,
  STAR = Star,
  TopHat2 = Tophat2
)

ColdataAll <- data.frame(
  sample    = rep(c("p15_s1","p15_s2","p15_s3","p5_s1","p5_s2","p5_s3"), 5),
  condition = rep(c("p15","p15","p15","p5","p5","p5"), 5),
  aligner   = rep(names(counts_list), each = 6)
)

CommonGene <- Reduce(intersect, map(counts_list, rownames))

# All PCA
CountAll <- counts_list |>
  map(~.x[CommonGene, ]) |>
  bind_cols() |>
  as.matrix()
colnames(CountAll) <- paste0(ColdataAll$aligner, "_", ColdataAll$sample)
rownames(ColdataAll) <- colnames(CountAll)
CountAll <- CountAll[rowSums(CountAll) > 15, ]

dds_all <- DESeqDataSetFromMatrix(
  countData = round(CountAll),
  colData   = ColdataAll,
  design    = ~ aligner + condition
)

dds_all <- dds_all[rowSums(counts(dds_all)) > 10, ]
vsd <- vst(dds_all, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = c("aligner", "condition"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

ggplot(pca_data, aes(x = PC1, y = PC2, color = aligner, shape = condition)) +
  geom_point(size = 5, alpha = 0.4) +
  scale_color_manual(values = c(
    BWA        = "tomato",
    HISAT2NOSS = "steelblue",
    HISAT2     = "lightgreen",
    STAR       = "gold",
    TopHat2    = "purple3"
  )) +
  scale_shape_manual(values = c(
    p5  = 16,  # 채운 원
    p15 = 15   # 채운 삼각형
  )) +
  labs(
    title = "PCA of All Samples across Aligners",
    x = paste0("PC1 (", pct_var[1], "%)"),
    y = paste0("PC2 (", pct_var[2], "%)")
  ) +
  theme_bw() +
  coord_fixed() +
  theme(aspect.ratio = 1)

# p5 PCA
counts_p5 <- CountAll[, ColdataAll$condition == "p5"]
coldata_p5 <- ColdataAll %>% filter(condition == "p5")

dds_p5 <- DESeqDataSetFromMatrix(
  countData = round(counts_p5),
  colData   = coldata_p5,
  design    = ~ aligner
)
dds_p5 <- dds_p5[rowSums(counts(dds_p5)) > 10, ]
vsd_p5 <- vst(dds_p5, blind = TRUE)

pca_p5  <- plotPCA(vsd_p5, intgroup = "aligner", returnData = TRUE)
pct_p5  <- round(100 * attr(pca_p5, "percentVar"), 1)

ggplot(pca_p5, aes(x = PC1, y = PC2, color = aligner)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_color_manual(values = c(
    BWA        = "tomato",
    HISAT2NOSS = "steelblue",
    HISAT2     = "lightgreen",
    STAR       = "gold",
    TopHat2    = "purple3"
  )) +
  labs(
    title = "PCA - p5 samples",
    x = paste0("PC1 (", pct_p5[1], "%)"),
    y = paste0("PC2 (", pct_p5[2], "%)")
  ) +
  theme_bw() +
  coord_fixed() +
  theme(aspect.ratio = 1)

# p15 PCA
counts_p15 <- CountAll[, ColdataAll$condition == "p15"]
coldata_p15 <- ColdataAll %>% filter(condition == "p15")

dds_p15 <- DESeqDataSetFromMatrix(
  countData = round(counts_p15),
  colData   = coldata_p15,
  design    = ~ aligner
)
dds_p15 <- dds_p15[rowSums(counts(dds_p15)) > 10, ]
vsd_p15 <- vst(dds_p15, blind = TRUE)

pca_p15  <- plotPCA(vsd_p15, intgroup = "aligner", returnData = TRUE)
pct_p15  <- round(100 * attr(pca_p15, "percentVar"), 1)

ggplot(pca_p15, aes(x = PC1, y = PC2, color = aligner)) +
  geom_point(size = 5, alpha = 0.8, shape = 15) +
  scale_color_manual(values = c(
    BWA        = "tomato",
    HISAT2NOSS = "steelblue",
    HISAT2     = "lightgreen",
    STAR       = "gold",
    TopHat2    = "purple3"
  )) +
  labs(
    title = "PCA - p15 samples",
    x = paste0("PC1 (", pct_p15[1], "%)"),
    y = paste0("PC2 (", pct_p15[2], "%)")
  ) +
  theme_bw() +
  coord_fixed() +
  theme(aspect.ratio = 1)
