options(iSEE.maxlevels=50)
stopifnot(suppressPackageStartupMessages({
    require(SingleCellExperiment)
    require(iSEE)
}))

# Load the (barely) preprocessed object ----

useHDF5 <- TRUE

if (useHDF5) {
    sce <- readRDS("sce.h5.rds")
} else {
    sce <- readRDS("sce.rds")
}

# Show the count matrix
assayNames(sce)
assay(sce, "counts")

# Show the gene metadata
rowData(sce)

# Show the cell metadata
colData(sce)

# Launch iSEE! ----

app <- iSEE(sce)

shiny::runApp(app, launch.browser = TRUE)
