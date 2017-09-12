
# Make sure we use scater prior to SingleCellExperiment release
stopifnot(packageVersion("scater") < package_version("1.5.11"))

# Load required packages ----

stopifnot(
  require(scater),
  require(Seurat),
  require(Matrix)
)

# Load data sets ----

# Normalised SCE with tSNE coordinates
sce.sc <- readRDS("rds/sce.sc.rds")
sce.sc

sce.mat <- Matrix(counts(sce.sc))

# Initialise Seurat ----

# Initialize the Seurat object with the raw (non-normalized data).
# Keep all genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected genes
scdc <- CreateSeuratObject(raw.data = sce.mat, min.cells = 3, min.genes = 200, project = "SC_DC")
scdc

# QC and filtering ----

# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
# For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
# We calculate the percentage of mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and non-log-normalized counts
# The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
# NOTE: You must have the Matrix package loaded to calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = scdc@data), value = TRUE)
percent.mito <- colSums(pbmc@raw.data[mito.genes, ]) / colSums(pbmc@raw.data)
