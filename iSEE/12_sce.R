
sce.endo_clusters <- readRDS("12_sce_infected/rds/sce.endo_clusters.rds")
sce.norm.tSNE <- readRDS("12_sce_infected/rds/sce.norm.tSNE.rds")

setdiff(
  colnames(colData(sce.endo_clusters)),
  colnames(colData(sce.norm.tSNE))
)

colData(sce.norm.tSNE) <- colData(sce.endo_clusters)[colnames(sce.endo_clusters),]

head(rownames(sce.norm.tSNE))

head(rowData(sce.norm.tSNE))

new_rownames <- with(
  rowData(sce.norm.tSNE),
  paste(gene_id, gene_name, sep = "__")
)

rownames(sce.norm.tSNE) <- new_rownames

saveRDS(sce.norm.tSNE, "iSEE/sce.rds")
