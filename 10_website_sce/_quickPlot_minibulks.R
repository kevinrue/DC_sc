
assayNames(sce.filtered)
assay(sce.filtered, "cpm")

geneIdFromGeneName <- function(x){
  subset(rowData(sce.filtered), gene_name == x, "gene_id", drop = TRUE)
}

# geneIdFromGeneName("MARCH1")

cpmFromGeneName <- function(x){
  geneId <- geneIdFromGeneName(x)
  cpmData <- log2(assay(sce.filtered, "cpm")[geneId,])
  res <- data.frame(
    cpm = cpmData,
    colData(sce.filtered)[,c("Time","Infection","Status","Donor")]
  )
  return(res)
}

# march1data <- cpmFromGeneName("MARCH1")

pointFromGeneName <- function(x){
  genedata <- cpmFromGeneName(x)
  ggplot(genedata) +
    facet_grid(Time ~ Donor) +
    geom_point(
      aes(interaction(Infection, Status), cpm, colour = Donor)
    ) +
    labs(
      title = x,
      x = NULL,
      y = expression(log[2]*" cpm")
    ) +
    theme(
      axis.text.x = element_text(angle = 90)
    )
}

lineFromGeneName <- function(x){
  genedata <- cpmFromGeneName(x)
  ggplot(genedata) +
    facet_grid(Donor ~ interaction(Infection, Status)) +
    geom_line(
      aes(Time, cpm, colour = interaction(Infection, Status), linetype = Donor, group = interaction(Infection, Status, Donor))
    ) +
    labs(
      title = x,
      x = NULL,
      y = expression(log[2]*" cpm")
    ) +
    theme(
      axis.text.x = element_text(angle = 90)
    ) +
    theme_bw()
}

pointFromGeneName("MARCH1")
lineFromGeneName("MARCH1")

pointFromGeneName("CTSL")
lineFromGeneName("CTSL")

pointFromGeneName("IL1B")

pointFromGeneName("IL12B")

pointFromGeneName("IL10")

pointFromGeneName("CTSS")

pointFromGeneName("CD83")

pointFromGeneName("CD86")

pointFromGeneName("HLA-DRB1")

pointFromGeneName("HLA-C")
