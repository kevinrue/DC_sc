# list.files("rds")

# Load the full object
sce <- readRDS("rds/sce.all.SCE.rds")

require(scater)
tmpSCE <- normalize(sce)
tmpSCE <- runPCA(tmpSCE)

# Fetch the reduced dimension data from the filtered & normalized object
tmpSCE <- readRDS("rds/sce.norm.tSNE.rds")
reducedDimNames(tmpSCE)
for (tmpName in reducedDimNames(tmpSCE)) {
    tmpSource <- reducedDim(tmpSCE, tmpName)
    tmpDest <- matrix(NA, nrow=ncol(sce), ncol=ncol(tmpSource), dimnames=list(colnames(sce), colnames(tmpSource)))
    tmpDest[colnames(tmpSCE), ] <- tmpSource
    reducedDim(sce, paste0("filtered_", tmpName)) <- tmpDest
}

library(iSEE)

# names(colData(sce))
colDataArgs <- colDataPlotDefaults(sce, 1)
colDataArgs$YAxis <- "Infection"
colDataArgs$XAxis <- "Column data"
colDataArgs$XAxisColData <- "Status"
colDataArgs$ColorBy <- "Column data"
colDataArgs$ColorByColData <- "Time"
colDataArgs$PointSize <- 2

app <- iSEE(sce, colDataArgs = colDataArgs)

library(shiny)
runApp(app)
