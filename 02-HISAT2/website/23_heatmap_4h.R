
library(scater)
library(scran)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
names(col.time) <- levels(sce.norm$Time)[1:3]
names(col.infection) <- levels(sce.norm$Infection)
names(col.status) <- levels(sce.norm$Status)

sce.endo <- sce.norm[!isSpike(sce.norm),]

# Colour map for heat map
exprsRange <- range(norm_exprs(sce.endo))
colorMap <- circlize::colorRamp2(
  c(min(exprsRange), median(exprsRange), max(exprsRange)),
  c("blue","white","red"))

# Subset by time point
sce.4h <- sce.endo[,sce.endo$Time == '4h']

# Subset genes
geneNames <- c("APOE", "APOC1", "LRP1", "CD1A", "CD1B", "CD1C", "CD1E")
geneNames <- c("APOE", "APOC1", "LRP1", "CD1A")
geneNames <- c("APOE", "CD1A")
geneMap <- with(
  fData(sce.4h),
  fData(sce.4h)[match(geneNames, gene_name), c("gene_name", "gene_id")]
  )
sce.genes <- sce.4h[geneMap$gene_id,]

# Refine color gradient
exprsRange <- range(norm_exprs(sce.genes)) * 1.25
colorMap <- circlize::colorRamp2(
  c(min(exprsRange), median(exprsRange), max(exprsRange)),
  c("blue","white","red"))

# Subset by Treatment
sce.mock <- sce.genes[,sce.genes$Infection == "Mock"]
mat.mock <- norm_exprs(sce.mock)

sce.d23inf <- sce.genes[,sce.genes$Treatment == "STM-D23580_Violet +"]
mat.d23inf <- norm_exprs(sce.d23inf)

sce.d23exp <- sce.genes[,sce.genes$Treatment == "STM-D23580_Exposed"]
mat.d23exp <- norm_exprs(sce.d23exp)

sce.lt2inf <- sce.genes[,sce.genes$Treatment == "STM-LT2_Violet +"]
mat.lt2inf <- norm_exprs(sce.lt2inf)

sce.lt2exp <- sce.genes[,sce.genes$Treatment == "STM-LT2_Exposed"]
mat.lt2exp <- norm_exprs(sce.lt2exp)
rownames(mat.lt2exp) <- geneNames

# pdf("23_out/heatmap_4h.pdf", height = 5, width = 8)
pdf("23_out/heatmap_4h_v2.pdf", height = 3, width = 8)
pdf("23_out/heatmap_4h_v3.pdf", height = 3, width = 8)
HeatmapList(
  Heatmap(
    matrix = mat.mock, col = colorMap, column_title = "Mock\nUninfected",
    show_row_names = FALSE, show_column_names = FALSE, cluster_rows = FALSE,
    show_heatmap_legend = FALSE) +
  Heatmap(
    matrix = mat.d23inf, col = colorMap, column_title = "STM-D23580\nViolet +",
    show_row_names = FALSE, show_column_names = FALSE,
    show_heatmap_legend = FALSE) +
  Heatmap(
    matrix = mat.d23exp, col = colorMap, column_title = "STM-D23580\nExposed",
    show_row_names = FALSE, show_column_names = FALSE,
    show_heatmap_legend = FALSE) +
  Heatmap(
    matrix = mat.lt2inf, col = colorMap, column_title = "STM-LT2\nViolet +",
    show_row_names = FALSE, show_column_names = FALSE,
    show_heatmap_legend = FALSE) +
  Heatmap(
    matrix = mat.lt2exp, col = colorMap, name = "norm.\nexpr.",
    column_title = "STM-LT2\nExposed",
    show_column_names = FALSE)
)
dev.off()
