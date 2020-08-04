
write.table(x = assay(sce, "counts"), file = "~/Desktop/DC/count_matrix.txt", quote = FALSE, sep = "\t")

write.table(colData(sce)[, c("Sample", "Infection", "Status", "Time", "Lane", "Plate",  "Well")], "~/Desktop/DC/cell_metadata.txt", quote = FALSE, sep = "\t")

write.table(rowData(sce)[, c("source", "type", "score", "phase", "gene_id", "gene_version",  "gene_name", "gene_source", "gene_biotype", "havana_gene", "havana_gene_version",  "transcript_id")], "~/Desktop/DC/gene_metadata.txt", quote = FALSE, sep = "\t")

