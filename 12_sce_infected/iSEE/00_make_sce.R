
stopifnot(suppressPackageStartupMessages({
    require(readxl)
    require(edgeR)
    require(rtracklayer)
    require(SingleCellExperiment)
    require(HDF5Array)
}))

# Import the cell metadata ----

# excel_sheets("samples.xlsx")
cellMetadataTable <- read_excel(path="samples.xlsx", sheet="pheno")

# Format various columns
cellMetadataTable$Infection <- factor(cellMetadataTable$Infection, c("Mock", "STM-LT2", "STM-D23580", "Blank"))
cellMetadataTable$Status <- factor(cellMetadataTable$Status, c("Uninfected", "Exposed", "Infected", "Bulk", "Blank"))
cellMetadataTable$Time <- factor(cellMetadataTable$Time, c("2h", "4h", "6h", "Blank"))
cellMetadataTable$Lane <- factor(cellMetadataTable$Lane)
cellMetadataTable$Plate <- factor(cellMetadataTable$Plate)
cellMetadataTable$Well <- factor(cellMetadataTable$Well)

# Add columns that combine other metadata
cellMetadataTable$Treatment <- with(cellMetadataTable,
    factor(x=droplevels(interaction(Infection, Status, sep="_"))))
cellMetadataTable$Group <- with(cellMetadataTable,
    factor(x=droplevels(interaction(Time, Treatment, sep="_"))))
cellMetadataTable$Sample <- with(cellMetadataTable,
    factor(x=droplevels(interaction(Group, Lane, Plate, Well, sep="_"))))

cellMetadataTable <- as.data.frame(cellMetadataTable, row.names=cellMetadataTable$Sample)

# Import the cell count data ----

RG <- with(cellMetadataTable,
    readDGE(files=File, path="counts", columns=c(1, 3), group=Group, labels=Sample))
# dim(RG)
# class(RG)

# Import gene metadata ----

gtfInfo <- import.gff("genes_with_ERCC.gtf")
names(gtfInfo) <- gtfInfo$gene_id

# Assemble the SingleCellExperiment object ----

sce <- SingleCellExperiment(
    assays=list(counts=RG$counts),
    colData=cellMetadataTable,
    rowRanges=gtfInfo)

h5file <- "sce.h5"
assay(sce, "counts") <- writeHDF5Array(assay(sce, "counts"), h5file, "counts", chunkdim = c(100, 100), verbose=TRUE)
saveRDS(file="sce.rds", sce)
