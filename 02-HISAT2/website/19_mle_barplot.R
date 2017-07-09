library(ggplot2)
library(dplyr)
library(RColorBrewer)

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")
scde.res <- readRDS("rds/scde.res_v8.rds")

# Colours ----

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
names(col.time) <- levels(sce.norm$Time)[1:3]
names(col.infection) <- levels(sce.norm$Infection)
names(col.status) <- levels(sce.norm$Status)

# Process ----

names(scde.res)

scde.x <- scde.res[["6h_STM-D23580_Violet +-6h_STM-LT2_Violet +"]]

head(scde.x)

scde.x$gene_name <- with(fData(sce.norm), gene_name[match(rownames(scde.x), gene_id)])

genenames.plot <- c(
  "STARD8", "BLOC1S3", "CLEC16A", "MARCH1", "USP1", "DNAJC13",
  "CLASP1", "TXNRD1", "TRIAP1", "CLECL1",
  "TFRC", "TRIM39", "VPS8", "USPL1", "RNF139", "LRR1", "NIT1", "NIT2"
)

data.frame(
  gene_name = genenames.plot,
  fData = genenames.plot %in% fData(sce.norm)[,"gene_name"]
)

scde.plot <- subset(scde.x, gene_name %in% genenames.plot)

scde.plot$gene_name <- with(scde.plot, reorder(gene_name, mle))

scde.plot$Highest <- with(scde.plot, ifelse(mle > 0, "STM-D23580", "STM-LT2"))

ggplot(scde.plot, aes(gene_name, mle)) +
  geom_col(aes(fill = Highest), alpha = 0.75) +
  coord_flip() +
  scale_fill_manual(values = col.infection) +
  scale_y_continuous(breaks = with(scde.plot, seq(floor(min(mle)), max(mle), 1))) +
  labs(
    x = NULL,
    y = "Maximum likelihood estimate of fold-change",
    fill = "Up-regulated in"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))
ggsave("19_out/mle_6h_D23-LT2_vertical.pdf", width = 7, height = 8)

ggplot(scde.plot, aes(gene_name, mle)) +
  geom_col(aes(fill = Highest), alpha = 0.75) +
  scale_fill_manual(values = col.infection) +
  scale_y_continuous(breaks = with(scde.plot, seq(floor(min(mle)), max(mle), 1))) +
  labs(
    x = NULL,
    y = "Maximum likelihood estimate of fold-change",
    fill = "Up-regulated in"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic"))
ggsave("19_out/mle_6h_D23-LT2_horizontal.pdf", width = 8, height = 6)
