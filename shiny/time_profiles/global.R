
stopifnot(
  requireNamespace("scater"),
  require(EnsDb.Hsapiens.v79),
  require(ggplot2),
  requireNamespace("dplyr")
)

workdir <- ifelse(interactive(), "../..", "../..")
folder.rds <- file.path(workdir, "rds")

sceset <- readRDS(file.path(folder.rds, "normalised.rds"))
featureNames(sceset) <- fData(sceset)$feature_id

annTable <- dplyr::arrange(
  ensembldb::select(
    EnsDb.Hsapiens.v79,
    fData(sceset)$feature_id,
    c("GENENAME", "GENEID"),
    "GENEID"
  ),
  GENENAME
)
