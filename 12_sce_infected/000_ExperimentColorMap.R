
library(iSEE)
library(RColorBrewer)

# col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
# names(col.status) <- levels(sce.norm$Status)
# sce.endo <- sce.norm[!isSpike(sce.norm),]
# sce.endo$TreatmentLabel <- with(colData(sce.endo), factor(
#   gsub("_", "\n", Treatment),
#   gsub("_", "\n", levels(Treatment))
# ))
# col.treatment <- RColorBrewer::brewer.pal(9, "Paired")[c(9,3,4,1,2)]
# names(col.treatment) <- gsub("_","\n",levels(sce.endo$TreatmentLabel))

timeColors <- function(n=3) {
    colors <- brewer.pal(9, "Set3")[c(2,7:8)]
    names(colors) <- c("2h", "4h", "6h")
    return(colors)
}

infectionColors <- function(n=3) {
    colors <- brewer.pal(12, "Paired")[c(9,4,2)]
    names(colors) <- c("Mock", "STM-LT2", "STM-D23580")
    return(colors)
}

quickClusterColor <- function(n=5L) {
    colors <- brewer.pal(8, "Dark2")[c(1, 2, 3, 4, 8)]
    names(colors) <- c(1, 2, 3, 0, "NA")
    return(colors)
}

ecm <- ExperimentColorMap(
    colData = list(
        Time=timeColors,
        Infection=infectionColors,
        quickCluster=quickClusterColor
    )
)
