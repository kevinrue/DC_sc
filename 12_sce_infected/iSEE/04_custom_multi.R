
stopifnot(suppressPackageStartupMessages({
    require(SingleCellExperiment)
    require(scran)
    require(scater)
    require(org.Hs.eg.db)
    require(RColorBrewer)
    require(iSEE)
}))

if (FALSE) {
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

    set.seed(1794)
    sce <- runTSNE(object = sce,
        ntop = 500,
        ncomponents = 2,
        pca = TRUE, initial_dims = 50,
        feature_set = which(rowData(sce)$source != "ERCC"),
        perplexity = round(mean(table(sce$Group)) / 2) # = 11
        )

    set.seed(1794)
    sce <- runDiffusionMap(object = sce,
        ntop = 500,
        ncomponents = 2,
        feature_set = which(rowData(sce)$source != "ERCC")
        )

    saveRDS(sce, "sce.04.rds")
}

sce <- readRDS("sce.04.rds")

# Preconfigure the initial state of the app ----

redDimArgs <- redDimPlotDefaults(sce, 1)

redDimArgs$Type <- c("TSNE")

redDimArgs$ColorBy <- "Column data"
redDimArgs$ColorByColData <- "Time"

redDimArgs$ShapeBy <- "Column data"
redDimArgs$ShapeByColData <- "Infection"

redDimArgs$PointSize <- 3

redDimArgs$LassoData[[1]] <- list(
    lasso = NULL, closed = TRUE,
    panelvar1 = NULL, panelvar2 = NULL,
    mapping = list(x = "X", y = "Y", colour = "ColorBy", shape = "ShapeBy"),
    coord = structure(c(0.96470530146804,
        1.07141748670796, 5.23319271106481, 8.32784608302247, 11.2090750845003, 9.07483137970191,
        2.67210026530675, 0.96470530146804, 20.4237156818852, 12.4695358942074, 11.2911388886255,
        15.1209291567667, 22.4859104416536, 24.5481052014219, 23.2224085701423, 20.4237156818852
        ), .Dim = c(8L, 2L)))

EMPTY_PLOT <- function() {
    ggplot() +
        theme_void() +
        geom_text(
            aes(x, y, label=label),
            data.frame(x=0, y=0, label="No column data selected."),
            size=5)
}

CUSTOM_PCA <- function(se, rows, columns, colour_by=NULL, size_by=NULL, scale_features=TRUE) {
    stopifnot(suppressMessages({require(scater)}))
    if (!is.null(columns)) {
        kept <- se[, columns]
    } else {
        return(EMPTY_PLOT())
    }

    scale_features <- as.logical(scale_features)
    kept <- runPCA(kept, feature_set=rows, scale_features=scale_features)
	plotPCA(kept, colour_by=colour_by, size_by = size_by)
}

CUSTOM_BARPLOT <- function(se, rows, columns) {
    if (!is.null(columns)) {
        kept <- se[, columns]
    } else {
        return(EMPTY_PLOT())
    }
    ggplot(as.data.frame(colData(kept))) +
        geom_bar(aes(gsub("_", "\n", Group), fill=gsub("_", "\n", Group))) +
        xlab("Group") +
        labs(fill="Group")
}

CUSTOM_MULTI <- function(se, rows, columns, colour_by=NULL, size_by=NULL, scale_features=TRUE) {
    stopifnot(suppressPackageStartupMessages(require(cowplot)))
    ggList <- list(
        PCA=CUSTOM_PCA(se, rows, columns, colour_by=colour_by, size_by=size_by, scale_features=scale_features),
        BARPLOT=CUSTOM_BARPLOT(se, rows, columns)
    )
    plot_grid(plotlist = ggList, ncol = 2)
}

customDataArgs <- customDataPlotDefaults(sce, 1)
customDataArgs$Function <- "CUSTOM_MULTI"
customDataArgs$Arguments <- "colour_by Treatment\nsize_by Infection\nscale_features TRUE"
customDataArgs$ColumnSource <- "Reduced dimension plot 1"
customDataArgs$SelectBoxOpen <- TRUE

initialPanels <- DataFrame(
    Name=c(
        "Reduced dimension plot 1",
        "Custom data plot 1"
    ),
    Width=c(4, 8)
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
    redDimArgs = redDimArgs, colDataArgs = colDataArgs, customDataArgs = customDataArgs,
    annotFun = annot.fun,
    initialPanels = initialPanels, colormap = ecm,
    customDataFun=list(CUSTOM_MULTI=CUSTOM_MULTI),
    appTitle = "Aulicino & Rue-Albrecht et al., 2018, Nat. Comm.")

# Launch iSEE! ----

shiny::runApp(app, launch.browser = TRUE)
