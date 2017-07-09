library(ggplot2)
library(RColorBrewer)

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
names(col.time) <- levels(sce.norm$Time)[1:3]
names(col.infection) <- levels(sce.norm$Infection)
names(col.status) <- levels(sce.norm$Status)

# MARCH1 ----

geneName <- "MARCH1"
geneId <- subset(fData(sce.norm), gene_name == geneName, gene_id, drop = TRUE)
cellsId <- subset(pData(sce.norm), Time == "6h", Sample, drop = TRUE)

gene.plot <- data.frame(
  norm_exprs = norm_exprs(sce.norm)[geneId, cellsId],
  pData(sce.norm)[cellsId, c("Time", "Infection", "Status")]
)

ggplot(
  gene.plot,
  aes(interaction(Infection, Status, sep = "\n"), norm_exprs)
  ) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  labs(
    title = sprintf("%s - %s", geneName, geneId),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20_out/violin_6h_%s.pdf", geneName), height = 4, width = 8)

# IL10 ----

geneName <- "IL10"
geneId <- subset(fData(sce.norm), gene_name == geneName, gene_id, drop = TRUE)
cellsId <- subset(pData(sce.norm), Time == "6h", Sample, drop = TRUE)

gene.plot <- data.frame(
  norm_exprs = norm_exprs(sce.norm)[geneId, cellsId],
  pData(sce.norm)[cellsId, c("Time", "Infection", "Status")]
)

ggplot(
  gene.plot,
  aes(interaction(Infection, Status, sep = "\n"), norm_exprs)
) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  labs(
    title = sprintf("%s - %s", geneName, geneId),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20_out/violin_6h_%s.pdf", geneName), height = 4, width = 8)
