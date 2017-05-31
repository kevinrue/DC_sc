sce.norm <- readRDS("data/sce.norm.tSNE.rds")

exprs_range <- range(norm_exprs(sce.norm))

geneNames <- sort(fData(sce.norm)[,"gene_name"])
geneNames <- geneNames[!is.na(geneNames)]

pairs_point <- function(
  data, mapping, ...,
  colour_map = NULL, colour_by = NULL, shape_by = NULL){

  p <- ggplot(data, mapping) +
    geom_point(aes_string(colour = colour_by, shape = shape_by), ...) +
    scale_x_continuous(limits = exprs_range) +
    scale_y_continuous(limits = exprs_range)

  if (!is.null(colour_map)){
    p <- p + scale_colour_manual(values = colour_map)
  }

  p
}

pairs_density <- function(
  data, mapping, ...,
  colour_map = NULL, fill_by = NULL){

  p <- ggplot(data, mapping) +
    geom_density(aes_string(fill = fill_by), ...) +
    scale_x_continuous(limits = exprs_range)

  if (!is.null(colour_map)){
    p <- p + scale_fill_manual(values = colour_map)
  }

  p
}
