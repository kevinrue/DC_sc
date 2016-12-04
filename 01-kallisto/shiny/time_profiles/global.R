
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

yMaxRange <- range(norm_exprs(sceset))

drawExpProfile <- function(
  exprMat, facetRow, facetCol, colourBy, shapeBy, yFullRange){

  g <- ggplot(exprMat, aes(x = Time, y = exprs))

  if (any(c(facetRow, facetCol) != ".")){
    f <- paste(facetRow, facetCol, sep = "~")
    message(f)
    g <- g + facet_grid(as.formula(f))
  }

  if (identical(colourBy, "None")){
    if (identical(shapeBy, "None")){
      g <- g +
        geom_violin(draw_quantiles = seq(0.25, 0.75, 0.25)) +
        geom_point()
    } else {
      g <- g +
        geom_violin(
          aes_string(linetype = shapeBy),
          draw_quantiles = seq(0.25, 0.75, 0.25)
        ) +
        geom_point(
          aes_string(shape = shapeBy),
          position = position_jitterdodge(0.25, 0)
        )
    }
  } else {
    if (identical(shapeBy, "None")){
      g <- g +
        geom_violin(
          aes_string(colour = colourBy),
          draw_quantiles = seq(0.25, 0.75, 0.25)
        ) +
        geom_point(
          aes_string(colour = colourBy),
          position = position_jitterdodge(0.25, 0, 0.9)
        )
    } else {
      g <- g +
        geom_violin(
          aes_string(colour = colourBy, linetype = shapeBy),
          draw_quantiles = seq(0.25, 0.75, 0.25)
        ) +
        geom_point(
          aes_string(colour = colourBy, shape = shapeBy),
          position = position_jitterdodge(0.25, 0, 0.9)
        )
    }
  }

  if (yFullRange){
    g <- g + scale_y_continuous(limits = yMaxRange)
  }

  g
}
