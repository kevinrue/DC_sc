
dir.create("32_out", showWarnings = FALSE)

# Make a simple SCE from the SCESet ----

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

rd_cols_idx <- which(
  !colnames(pData(sce.norm@featureData)) %in%
    c("seqnames","ranges", "strand", "start", "end", "width", "element"))
gr <- GRanges(DataFrame(pData(sce.norm@featureData))[,-rd_cols_idx])
mcols(gr) <- DataFrame(pData(sce.norm@featureData))[,rd_cols_idx]

sce.norm.SCE <- SingleCellExperiment(
  assays = SimpleList(
    counts = sce.norm@assayData$counts,
    norm_exprs = sce.norm@assayData$norm_exprs
    ),
  rowRanges = gr,
  colData = DataFrame(pData(sce.norm@phenoData))
)

saveRDS(sce.norm.SCE, "rds/sce.norm.SCE_32.rds")

# Error in validObject(.Object) :
#   invalid class “GRangesList” object: 'mcols(x)' cannot have columns named "seqnames",  "ranges", "strand", "start", "end", "width",  or "element"

# Flag ERCC spike-in features ----

erccIds <- grep("^ERCC", rownames(sce.norm.SCE), value = TRUE)
isSpike(sce.norm.SCE, "ERCC") <- erccIds
table(isSpike(sce.norm.SCE))

# Remove spike-in features ----

sce.endo <- sce.norm.SCE[!isSpike(sce.norm.SCE),]
dim(sce.endo)

# Prepare colour maps for the various phenotypes ----

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
names(col.time) <- levels(sce.endo$Time)[1:3]
names(col.infection) <- levels(sce.endo$Infection)
names(col.status) <- levels(sce.endo$Status)

# Prepare colour map for heat map ----

exprsRange_all <- range(assay(sce.endo, "norm_exprs"))
colorMap <- circlize::colorRamp2(
  c(min(exprsRange_all), median(exprsRange_all), max(exprsRange_all)),
  c("blue","white","red"))

# Map gene names to IDs ----

mapGeneNames <- function(query, x){
  stopifnot(length(query) > 0, is.character(query))
  stopifnot(!any(duplicated(query)))
  # res <- subset(
  #   x=rowData(x),
  #   subset=gene_name %in% query,
  #   select=c("gene_name","gene_id"))
  res <- with(
    rowData(x),
    rowData(x)[match(query, gene_name), c("gene_name", "gene_id")]
  )
  dupGeneNames <- duplicated(res$gene_name)
  if (any(dupGeneNames)){
    warning("Gene names matched to multiple IDs:")
    print(res[duplicated(res$gene_name),])
  }
  missedGenes <- setdiff(query, res$gene_name)
  if (length(missedGenes) > 0){
    stop(
      paste(c("Gene names not found:", missedGenes), collapse = "\n"))
  }
  return(res)
}

# genes_2h <- c("DHRS9", "PLA2G16", "ZDHHC8", "GOT1")
# mapGeneNames(genes_2h, sce.endo)

# Subset time point ----

subsetTime <- function(query, x){
  query <- match.arg(query, levels(x$Time), several.ok = FALSE)
  return(x[,x$Time == query])
}

# subsetTime("4h", sce.endo)$Time

# Heat map ----

drawHeatMap <- function(geneNames, x, time){
  # Subset to genes
  sce.genes <- x[mapGeneNames(geneNames, x)$gene_id,]
  stopifnot(nrow(sce.genes) > 0)

  # Refine colour map for expression
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

  return(
    HeatmapList(
      Heatmap(
        matrix = mat.mock, col = colorMap,
        column_title = "Mock\nUninfected",
        show_row_names = FALSE, show_column_names = FALSE,
        cluster_rows = FALSE,
        show_heatmap_legend = FALSE) +
        Heatmap(
          matrix = mat.d23inf, col = colorMap,
          column_title = "STM-D23580\nViolet +",
          show_row_names = FALSE, show_column_names = FALSE,
          show_heatmap_legend = FALSE) +
        Heatmap(
          matrix = mat.d23exp, col = colorMap,
          column_title = "STM-D23580\nExposed",
          show_row_names = FALSE, show_column_names = FALSE,
          show_heatmap_legend = FALSE) +
        Heatmap(
          matrix = mat.lt2inf, col = colorMap,
          column_title = "STM-LT2\nViolet +",
          show_row_names = FALSE, show_column_names = FALSE,
          show_heatmap_legend = FALSE) +
        Heatmap(
          matrix = mat.lt2exp, col = colorMap,
          name = "norm.\nexpr.",
          column_title = "STM-LT2\nExposed",
          show_column_names = FALSE)
    )
  )
}

drawHeatMap(c("DHRS9","PLA2G16","ZDHHC8","GOT1"), sce.endo, "2h")
pdf("32_out/2h_01_DHRS9.pdf", height = 2, width = 8)
drawHeatMap(c("DHRS9"), sce.endo, "2h")
dev.off()

drawHeatMap(
  c("ELOVL7","FABP5","APOL1","APOL2","APOL6","ACSL3","HCAR2","ABHD6",
    "CH25H","NR1H3","FFAR4"),
  sce.endo, "4h")

pdf("32_out/4h_01.pdf", height = 5, width = 8)
drawHeatMap(
  c("ELOVL7","FABP5","APOL1","APOL2","APOL6","ACSL3","HCAR2","ABHD6",
    "NR1H3"),
  sce.endo, "4h")
dev.off()
