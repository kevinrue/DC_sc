library(ggplot2)
library(RColorBrewer)

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

col.time <- brewer.pal(9, "Set3")[c(2,7:8)]
col.infection <- brewer.pal(12, "Paired")[c(9,2,4)]
col.status <- brewer.pal(12, "Paired")[c(9,1,7)]
names(col.time) <- levels(sce.norm$Time)[1:3]
names(col.infection) <- levels(sce.norm$Infection)
names(col.status) <- levels(sce.norm$Status)

exprsRange <- range(norm_exprs(sce.norm))

# APOE ----

geneName <- "APOE"
geneId <- subset(fData(sce.norm), gene_name == geneName, gene_id, drop = TRUE)

gene.plot <- data.frame(
  norm_exprs = norm_exprs(sce.norm)[geneId, ],
  pData(sce.norm)[, c("Time", "Infection", "Status")]
)

ggplot(
  gene.plot,
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  facet_grid(Time ~ .) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s", geneName, geneId),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20_out/violin_all_%s.pdf", geneName), height = 6, width = 8)

t_h <- "4h"

ggplot(
  subset(gene.plot, Time == t_h),
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s (%s)", geneName, geneId, t_h),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20_out/violin_%s_%s.pdf", t_h, geneName), height = 4, width = 8)


# PLAT ----

geneName <- "PLAT"
geneId <- subset(fData(sce.norm), gene_name == geneName, gene_id, drop = TRUE)

gene.plot <- data.frame(
  norm_exprs = norm_exprs(sce.norm)[geneId, ],
  pData(sce.norm)[, c("Time", "Infection", "Status")]
)

ggplot(
  gene.plot,
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  facet_grid(Time ~ .) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s", geneName, geneId),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20b_out/violin_all_%s.pdf", geneName), height = 6, width = 8)

t_h <- "4h"

ggplot(
  subset(gene.plot, Time == t_h),
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s (%s)", geneName, geneId, t_h),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20b_out/violin_%s_%s.pdf", t_h, geneName), height = 4, width = 8)

# EBI3 ----

geneName <- "EBI3"
geneId <- subset(fData(sce.norm), gene_name == geneName, gene_id, drop = TRUE)

gene.plot <- data.frame(
  norm_exprs = norm_exprs(sce.norm)[geneId, ],
  pData(sce.norm)[, c("Time", "Infection", "Status")]
)

ggplot(
  gene.plot,
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  facet_grid(Time ~ .) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s", geneName, geneId),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20b_out/violin_all_%s.pdf", geneName), height = 6, width = 8)

t_h <- "4h"

ggplot(
  subset(gene.plot, Time == t_h),
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s (%s)", geneName, geneId, t_h),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20b_out/violin_%s_%s.pdf", t_h, geneName), height = 4, width = 8)

# CTSL ----

geneName <- "CTSL"
geneId <- subset(fData(sce.norm), gene_name == geneName, gene_id, drop = TRUE)

gene.plot <- data.frame(
  norm_exprs = norm_exprs(sce.norm)[geneId, ],
  pData(sce.norm)[, c("Time", "Infection", "Status")]
)

ggplot(
  gene.plot,
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  facet_grid(Time ~ .) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s", geneName, geneId),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20b_out/violin_all_%s.pdf", geneName), height = 6, width = 8)

t_h <- "4h"

ggplot(
  subset(gene.plot, Time == t_h),
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s (%s)", geneName, geneId, t_h),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20b_out/violin_%s_%s.pdf", t_h, geneName), height = 4, width = 8)

# MS4A4A ----

geneName <- "MS4A4A"
geneId <- subset(fData(sce.norm), gene_name == geneName, gene_id, drop = TRUE)

gene.plot <- data.frame(
  norm_exprs = norm_exprs(sce.norm)[geneId, ],
  pData(sce.norm)[, c("Time", "Infection", "Status")]
)

ggplot(
  gene.plot,
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  facet_grid(Time ~ .) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s", geneName, geneId),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20b_out/violin_all_%s.pdf", geneName), height = 6, width = 8)

t_h <- "4h"

ggplot(
  subset(gene.plot, Time == t_h),
  aes(interaction(Status, Infection, sep = "\n"), norm_exprs)
) +
  geom_violin(
    aes(fill = Infection),
    draw_quantiles = c(0.25, 0.5, 0.75),
    alpha = 0.5
  ) +
  geom_jitter(aes(colour = Infection), width = 0.3, height = 0) +
  scale_fill_manual(values = col.infection) +
  scale_colour_manual(values = col.infection) +
  scale_y_continuous(limits = exprsRange) +
  labs(
    title = sprintf("%s - %s (%s)", geneName, geneId, t_h),
    x = NULL,
    y = "Normalised expression"
  ) +
  theme_minimal()

ggsave(sprintf("20b_out/violin_%s_%s.pdf", t_h, geneName), height = 4, width = 8)
