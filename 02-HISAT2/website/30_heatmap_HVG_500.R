
# Make sure we use scater prior to SingleCellExperiment release
stopifnot(packageVersion("scater") < package_version("1.5.11"))

# Load required packages ----

stopifnot(
  require(ComplexHeatmap),
  require(RColorBrewer),
  requireNamespace("circlize")
)

# Load data sets ----

# Normalised SCE with tSNE coordinates
sce.norm <- readRDS("rds/sce.norm.tSNE.rds")
sce.norm

# Table of HVGs
hvg.out <- read.csv("07_out/hvg.out.csv", row.names = 1)
head(hvg.out)

# Prepare colours ----

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
names(col.time) <- levels(sce.norm$Time)
names(col.infection) <- levels(sce.norm$Infection)
names(col.status) <- levels(sce.norm$Status)

# Heat map expression levels
exprsRange <- range(norm_exprs(sce.norm))
colorMap <- circlize::colorRamp2(
  c(min(exprsRange), median(exprsRange)*0.9, max(exprsRange)*0.9),
  c("black","purple","yellow")
  # c("blue","white","red")
)

# Preprocess data set ----

# Order by decreasing biological variance
hvg.out <- hvg.out[order(hvg.out$bio, decreasing = TRUE),]
sce.hvg <- sce.norm[rownames(hvg.out),]
sce.hvg

# There shouldn't be any ERCC left, but filter anyway
isSpike(sce.hvg)
table(isSpike(sce.hvg))
sce.endo <- sce.hvg[!isSpike(sce.hvg),]
sce.endo

# Take top N features (by decreasing bio variance)
sce.top <- sce.endo[1:500,]
sce.top

# Make heat map ----

# Cluster HVGs by expression profile over all cells
hvg.endo.d <- dist(norm_exprs(sce.top))
hvg.endo.clust <- hclust(hvg.endo.d)
rm(hvg.endo.d)

# Cluster samples by expression within each experimental group
sample.order <- c()
expGroups <- unique(pData(sce.top)[,c("Time","Infection","Status")])
expGroups <- dplyr::arrange(expGroups, Time, Status, Infection)
for (groupIndex in seq_len(nrow(expGroups))){
  time <- expGroups$Time[groupIndex]
  infection <- expGroups$Infection[groupIndex]
  status <- expGroups$Status[groupIndex]
  sample.index <- with(pData(sce.top), which(
    Time == time & Infection == infection & Status == status
  ))
  hvg.sample.d <- dist(t(norm_exprs(sce.top[,sample.index])))
  local.order <- sample.index[hclust(hvg.sample.d)$order]
  sample.order <- c(sample.order, local.order)
  rm(sample.index, local.order, time, infection, status)
}
rm(expGroups, groupIndex)

sce.reorder <- sce.top[,sample.order]

# Split samples by time point
sce.2h <- sce.reorder[,sce.reorder$Time == "2h"]
sce.4h <- sce.reorder[,sce.reorder$Time == "4h"]
sce.6h <- sce.reorder[,sce.reorder$Time == "6h"]

# Draw heat map annotations from phenotype information
ha_2h <- HeatmapAnnotation(
  df = pData(sce.2h)[,c("Time", "Status", "Infection")],
  col = list(
    Time = col.time,
    Status = col.status,
    Infection = col.infection
  )
)
ha_4h <- HeatmapAnnotation(
  df = pData(sce.4h)[,c("Time", "Status", "Infection")],
  col = list(
    Time = col.time,
    Status = col.status,
    Infection = col.infection
  ), show_legend = FALSE
)
ha_6h <- HeatmapAnnotation(
  df = pData(sce.6h)[,c("Time", "Status", "Infection")],
  col = list(
    Time = col.time,
    Status = col.status,
    Infection = col.infection
  ), show_legend = FALSE
)

# Draw main heat maps
ht_2h <- Heatmap(
  norm_exprs(sce.2h), col = colorMap,
  name = "norm.\nexprs.", column_title = "2h",
  top_annotation = ha_2h,
  row_order = hvg.endo.clust$order, column_order = NULL,
  cluster_rows = hvg.endo.clust, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE, width = unit(4, "inches")
)
ht_4h <- Heatmap(
  norm_exprs(sce.4h), col = colorMap,
  name = "4h", column_title = "4h",
  top_annotation = ha_4h,
  row_order = hvg.endo.clust$order, column_order = NULL,
  cluster_rows = hvg.endo.clust, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  show_heatmap_legend = FALSE, width = unit(4, "inches")
)
ht_6h <- Heatmap(
  norm_exprs(sce.6h), col = colorMap,
  name = "6h", column_title = "6h",
  top_annotation = ha_6h,
  row_order = hvg.endo.clust$order, column_order = NULL,
  cluster_rows = hvg.endo.clust, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  show_heatmap_legend = FALSE , width = unit(4, "inches")
)

# Draw heat map ----

pdf("30_out/heatmap_split_time.pdf", height = 6, width = 12)
draw(ht_2h + ht_4h + ht_6h)
dev.off()

# Cut the tree of 500 into blocks ----

hvg.endo.clust
class(hvg.endo.clust)
plot(hvg.endo.clust)
abline(h=100, col="red")

h_100 <- cutree(hvg.endo.clust, h = 100)
table(h_100) # 33 clusters at height 100

k_10 <- cutree(hvg.endo.clust, k = 10)
table(k_10)

# Draw heat map augmented with cluster membership ----

ht_tree <- Heatmap(
  sprintf("%2d", k_10), col = brewer.pal(10, "Paired"),
  name = "k", column_title = "k",
  row_order = hvg.endo.clust$order, column_order = NULL,
  cluster_rows = hvg.endo.clust, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  width = unit(1/2, "inches")
)
pdf("30_out/heatmap_split_time_k.pdf", height = 6, width = 15)
draw(ht_2h + ht_4h + ht_6h + ht_tree)
dev.off()

# Select representativ genes by cluster ----

# Table of cluster membership for each gene
gene_k <- data.frame(
  fData(sce.6h)[
    hvg.endo.clust$order,c("gene_id","gene_name")],
  k_10 = k_10[hvg.endo.clust$order]
)
write.csv(gene_k, "30_out/gene_k_10.csv", row.names = FALSE)

# Annotate selected genes from each cluster -----

# Test: genes from cluster 4
k_10.4 <- c("IRF8","OAS2","OAS3")

# anno.names <- c(
#   k_10.4
# )

anno.xlsx <- "30_out/gene_k_10.xlsx"

anno.data <- xlsx::read.xlsx(anno.xlsx, sheetName = "gene_k_10")
dim(anno.data)

anno.subset <- subset(anno.data, Heatmap == "x", select=c("gene_name","gene_id"))
dim(anno.subset)

gene.idx = match(anno.subset$gene_id, rownames(sce.6h))

ra <- rowAnnotation(
  link = row_anno_link(
    at = gene.idx, labels = as.character(anno.subset$gene_name),
    labels_gp = gpar(cex=2/3), link_width = unit(0.5, "inches")
    ), width = max_text_width(labels)
)

pdf("30_out/heatmap_split_time_anno2.pdf", height = 10, width = 17)
# draw(ht_2h + ht_4h + ht_6h + ra)
draw(ht_2h + ht_4h + ht_6h + ht_tree + ra)
dev.off()

