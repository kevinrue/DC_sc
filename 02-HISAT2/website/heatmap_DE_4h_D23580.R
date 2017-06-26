
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(broom)

# Normalised SCESet
sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

# Endogenous genes
sce.endo <- sce.norm[!isSpike(sce.norm),]

# Colour range for heat map
exprsRange <- range(norm_exprs(sce.endo))
colorMap.exprs <- colorRamp2(
  c(min(exprsRange), median(exprsRange), max(exprsRange)),
  c("blue","white","red"))

# Colour code for phenotypes
col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
names(col.time) <- levels(sce.endo$Time)

col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
names(col.infection) <- levels(sce.endo$Infection)

col.status <- brewer.pal(12, "Paired")[c(9,7,1)]
names(col.status) <- levels(sce.endo$Status)

col.cluster <- brewer.pal(9, "Set3")[c(3,4,5,6,9)]
names(col.cluster) <- c(1,2,3,0,"NA")

# Normalised expression matrix
# ne <- norm_exprs(sce.endo)

# Selected gene names
geneNames <- c(
  "APOE", "APOC1", "PI4KB", "HADHA", "LRP1", "NR1H3", "OLR1",
  "CST3", "HLA-DQB1", "HLA-DQA2", "CCL23", "TBC1D17", "TRIM21", "CLEC4G"
)
ht.split <- rep(c(1, 2), c(7, 7))

# Matching gene identifiers
geneIds <- with(fData(sce.endo), gene_id[match(geneNames, gene_name)])

exprs.d23inf <- norm_exprs(sce.endo)[
  geneIds, sce.endo$Time == "4h" & sce.endo$Treatment == "STM-D23580_Violet +"
]
exprs.d23exp <- norm_exprs(sce.endo)[
  geneIds, sce.endo$Time == "4h" & sce.endo$Treatment == "STM-D23580_Exposed"
  ]
exprs.lt2inf <- norm_exprs(sce.endo)[
  geneIds, sce.endo$Time == "4h" & sce.endo$Treatment == "STM-LT2_Violet +"
  ]
exprs.lt2exp <- norm_exprs(sce.endo)[
  geneIds, sce.endo$Time == "4h" & sce.endo$Treatment == "STM-LT2_Exposed"
  ]

rownames(exprs.d23exp) <- geneNames
rownames(exprs.d23inf) <- geneNames
rownames(exprs.lt2exp) <- geneNames
rownames(exprs.lt2inf) <- geneNames

ht_list <- Heatmap(
  exprs.d23inf, colorMap.exprs, name = "n.e.", column_title = "STM-D23580 Violet +",
  top_annotation = new("HeatmapAnnotation"), row_names_side = "left",
  cluster_rows = FALSE, cluster_columns = TRUE,
  show_row_dend = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
  split = ht.split
  ) +
  Heatmap(
    exprs.d23exp, colorMap.exprs, name = "n.e.2", column_title = "STM-D23580 Exposed",
    top_annotation = new("HeatmapAnnotation"), row_names_side = "left",
    cluster_rows = FALSE, cluster_columns = TRUE,
    show_row_dend = FALSE, show_column_dend = TRUE,
    show_row_names = FALSE, show_column_names = FALSE,
    show_heatmap_legend = FALSE
  ) +
  Heatmap(
    exprs.lt2inf, colorMap.exprs, name = "n.e.3", column_title = "STM-LT2 Violet +",
    top_annotation = new("HeatmapAnnotation"), row_names_side = "left",
    cluster_rows = FALSE, cluster_columns = TRUE,
    show_row_dend = FALSE, show_column_dend = TRUE,
    show_row_names = FALSE, show_column_names = FALSE,
    show_heatmap_legend = FALSE
  ) +
  Heatmap(
    exprs.lt2exp, colorMap.exprs, name = "n.e.4", column_title = "STM-LT2 Exposed",
    top_annotation = new("HeatmapAnnotation"), row_names_side = "left",
    cluster_rows = FALSE, cluster_columns = TRUE,
    show_row_dend = FALSE, show_column_dend = TRUE,
    show_row_names = FALSE, show_column_names = FALSE,
    show_heatmap_legend = FALSE
  )

pdf(file = "heatmap_DE_4h_D23580_out/APOE_correlated.pdf", width = 10, height = 4)
draw(ht_list, gap = unit(10, "mm"))
dev.off()

# Gene inversely correlated with APOE ----

APOE.id <- with(fData(sce.endo), gene_id[match("APOE", gene_name)])

APOE.cor <- do.call("rbind", apply(
  norm_exprs(sce.endo)[
    which(featureNames(sce.endo) != APOE.id),
    sce.endo$Time == "4h" & sce.endo$Treatment == "STM-D23580_Violet +"
  ],
  1,
  function(x){
    tidy(cor.test(
      norm_exprs(sce.endo)[
        APOE.id,
        sce.endo$Time == "4h" & sce.endo$Treatment == "STM-D23580_Violet +"
      ],
      x,
      method = "pearson"))
  }
))

APOE.cor$BH <- p.adjust(APOE.cor$p.value, "BH")
APOE.cor$gene_name <-
  with(fData(sce.endo), gene_name[match(rownames(APOE.cor), gene_id)])
