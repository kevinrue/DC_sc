
sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

col.time <- brewer.pal(3, "Accent")[c(1,2,3)]
col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
names(col.time) <- levels(sce.norm$Time)

col.infection <- brewer.pal(5, "Set1")[c(5,2,3)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
names(col.infection) <- levels(sce.norm$Infection)

col.status <- brewer.pal(8, "Set1")[c(5,8,4)]
col.status <- brewer.pal(12, "Paired")[c(9,7,1)]
names(col.status) <- levels(sce.norm$Status)

col.cluster <- brewer.pal(9, "Set3")[c(3,4,5,6,9)]
names(col.cluster) <- c(1,2,3,0,"NA")

# Heat map marker genes ----

exprsRange <- range(norm_exprs(sce.norm))
colorMap.exprs <- colorRamp2(
  c(min(exprsRange), median(exprsRange), max(exprsRange)),
  c("blue","white","red"))

sce.2h <- sce.norm[,sce.norm$Time == "2h"]
sce.4h <- sce.norm[,sce.norm$Time == "4h"]
sce.6h <- sce.norm[,sce.norm$Time == "6h"]

d.e.2h <- dist(reducedDimension(sce.2h), diag = TRUE, upper = TRUE)
mat.e.2h <- as.matrix(d.e.2h)
h.e.2h <- hclust(d.e.2h)
ord.mat.2h <- mat.e.2h[h.e.2h$order, h.e.2h$order]
rm(d.e.2h, ord.mat.2h)

geneNames.ht <- list(
  "1. Trans. reprog."=c("DUSP1","DUSP2","SOCS3","IER3","NFKBIZ"),
  "2. Interm."=c("IL1B","IL1A","IFNB1","CCL7","TNF"),
  "3. Late"=c("IL6","IL10","IL12B","IL27")
)
geneIds.ht <- with(
  fData(sce.norm),
  gene_id[match(unlist(geneNames.ht), gene_name)]
)

ht_column <- HeatmapAnnotation(
  df = pData(sce.2h)[,c("Time", "Status", "Infection")],
  col = list(Time = col.time, Status = col.status, Infection = col.infection)
)
exprs.2h <- norm_exprs(sce.2h)[geneIds.ht,]
rownames(exprs.2h) <- unlist(geneNames.ht)
ht_gene_2h <- Heatmap(
  exprs.2h, colorMap.exprs, name = "2h", column_title = "2h",
  top_annotation = ht_column, row_names_side = "left",
  cluster_rows = FALSE, cluster_columns = h.e.2h,
  show_row_dend = FALSE, show_column_names = FALSE,
  split = rep(names(geneNames.ht), times = lengths(geneNames.ht))
); ht_gene_2h
