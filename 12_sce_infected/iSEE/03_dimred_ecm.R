
stopifnot(suppressPackageStartupMessages({
    require(SingleCellExperiment)
    require(scran)
    require(scater)
    require(org.Hs.eg.db)
    require(RColorBrewer)
    require(iSEE)
}))

# Load the (barely) preprocessed object ----

sce <- readRDS("sce.rds")

# Filter to the good cells ----

goodCellNames <- scan(file = "good_cells.txt", what = "character")
sce <- sce[, goodCellNames]

# Clean up empty factor levels
colData(sce) <- droplevels(colData(sce))

# Normalize ----

# Compute broad cluster membership
sce$quickCluster <- quickCluster(sce, min.size=min(table(sce$Group)))
table(sce$quickCluster)

# Compute sum factors using cell pools within broad cluster
sce <- computeSumFactors(sce, sizes = c(10, 15), clusters = sce$quickCluster)

# Normalize using the size factors computed above
sce <- normalize(sce)

# Compute dimensionality reduction results ----

sce <- runPCA(object = sce,
    ntop = 500,
    ncomponents = 10,
    feature_set = which(rowData(sce)$source != "ERCC"))

sce <- runTSNE(object = sce,
    ntop = 500,
    ncomponents = 2,
    pca = TRUE, initial_dims = 50,
    feature_set = which(rowData(sce)$source != "ERCC"),
    perplexity = round(mean(table(sce$Group)) / 2) # = 11
    )

sce <- runDiffusionMap(object = sce,
    ntop = 500,
    ncomponents = 2,
    feature_set = which(rowData(sce)$source != "ERCC")
    )

# Preconfigure the initial state of the app ----

redDimArgs <- redDimPlotDefaults(sce, 3)

redDimArgs$Type <- c("PCA", "TSNE", "DiffusionMap")

redDimArgs$ColorBy <- "Column data"
redDimArgs$ColorByColData <- "Time"

redDimArgs$ShapeBy <- "Column data"
redDimArgs$ShapeByColData <- "Status"

redDimArgs$PointSize <- 3

redDimArgs$SelectByPlot <- "Column data plot 1"
redDimArgs$SelectAlpha <- 0.3

colDataArgs <- colDataPlotDefaults(sce, 1)

colDataArgs$YAxis <- "Infection"

colDataArgs$XAxis <- "Column data"
colDataArgs$XAxisColData <- "Status"

colDataArgs$ColorBy <- "Column data"
colDataArgs$ColorByColData <- "Infection"

redDimArgs$PointSize <- 3

colDataArgs$BrushData[[1]] <- list(
    xmin = 0.5, xmax = 1.5, ymin = 0.5, ymax = 1.5, mapping = list(x = "X", y = "Y"),
    log = list(x = NULL, y = NULL), direction = "xy",
    brushId = "colDataPlot1_Brush", outputId = "colDataPlot1")

initialPanels <- DataFrame(
    Name=c(
        "Reduced dimension plot 1",
        "Reduced dimension plot 2",
        "Reduced dimension plot 3",
        "Column data plot 1"
    ),
    Width=c(4, 4, 4, 6)
)

# gene_biotype has 46 unique values (including NA)
# however, it should be treated as a categorical covariate
options(iSEE.maxlevels=50)

# Setting up the annotation function ----

annot.fun <- annotateEnsembl(sce, orgdb=org.Hs.eg.db, keytype="ENSEMBL", rowdata_col="gene_id")

# Setting up colormaps ----

colormapInfection <- function(n) {
    x <- brewer.pal(12, "Paired")[c(9, 4, 2, 5)]
    names(x) <- c("Mock", "STM-LT2", "STM-D23580", "Blank")
    x
}

colormapStatus <- function(n) {
    x <- brewer.pal(12, "Paired")[c(9, 1, 7, 3, 5)]
    names(x) <- c("Uninfected", "Exposed", "Infected", "Bulk", "Blank")
    x
}

colormapTime <- function(n) {
    x <- brewer.pal(9, "Set3")[c(2, 7, 8, 9)]
    names(x) <- c("2h", "4h", "6h", "Blank")
    x
}

ecm <- ExperimentColorMap(
    colData=list(
        Infection=colormapInfection,
        Status=colormapStatus,
        Time=colormapTime
    )
)

app <- iSEE(
    se = sce,
    redDimArgs = redDimArgs, colDataArgs = colDataArgs,
    annotFun = annot.fun,
    initialPanels = initialPanels, colormap = ecm,
    appTitle = "Aulicino & Rue-Albrecht et al., 2018, Nat. Comm.")

# Launch iSEE! ----

shiny::runApp(app, launch.browser = TRUE)
