sce <- readRDS("12_sce_infected/rds/sce.norm.rds")

write.table(x = assay(sce, "counts"), file = "12_sce_infected/export/count_matrix.txt", quote = FALSE, sep = "\t")
write.table(x = assay(sce, "logcounts"), file = "12_sce_infected/export/logcount_matrix.txt", quote = FALSE, sep = "\t")

write.table(colData(sce)[, c("Sample", "Infection", "Status", "Time", "Lane", "Plate",  "Well")], "12_sce_infected/export/cell_metadata.txt", quote = FALSE, sep = "\t")

write.table(rowData(sce)[, c("source", "type", "score", "phase", "gene_id", "gene_version",  "gene_name", "gene_source", "gene_biotype", "havana_gene", "havana_gene_version",  "transcript_id")], "12_sce_infected/export/gene_metadata.txt", quote = FALSE, sep = "\t")

