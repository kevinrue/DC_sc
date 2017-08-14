
library(ggplot2)
library(scater)

# Load data ----

sce.norm <- readRDS("rds/sce.norm.tSNE.rds")
sce.endo <- sce.norm[!isSpike(sce.norm),]

# IL1B feature identifier ----

il1b.id <- with(fData(sce.endo), gene_id[gene_name == "IL1B"])

il1b.data <- data.frame(
  IL1B = norm_exprs(sce.endo)[il1b.id,],
  tSNE1 = reducedDimension(sce.endo)[,1],
  tSNE2 = reducedDimension(sce.endo)[,2],
  pData(sce.endo)[,c("Time","Infection","Status","Treatment")]
)

# RColorBrewer::display.brewer.all()
# col.reds <- RColorBrewer::brewer.pal(6, "Paired")[5:6]
col.gradient <- RColorBrewer::brewer.pal(11, "RdYlBu")
# col.br <- RColorBrewer::brewer.pal(2, "Set1")[2:1]
col.treatment <- RColorBrewer::brewer.pal(9, "Paired")[c(1:4,9)]
names(col.treatment) <- c(paste(
  rep(
    paste("STM", c("D23580", "LT2"), sep = "-"),
    each = 2
  ),
  c("Exposed", "Violet +")
), "Mock Uninfected")

ggplot(il1b.data, aes(tSNE1, tSNE2, shape = paste(Infection, Status))) +
  geom_rug() +
  geom_point(aes(colour = IL1B), alpha = 0.75, size = 1.5) +
  stat_ellipse(
    aes(
      fill = paste(Infection, Status)
    ),
    subset(il1b.data, Time == "6h"),
    alpha = 0.4,
    level = 0.5,
    geom = "polygon"
    # geom = "path"
  ) +
  scale_fill_manual(values = col.treatment) +
  labs(fill = "Treatment", shape = "Treatment") +
  theme_bw()
ggsave("28_out/tSNE_ellipses_fill.pdf", height = 5, width = 8)

# ggplot(il1b.data, aes(tSNE1, tSNE2, shape = paste(Infection, Status))) +
#   geom_rug() +
#   geom_point(aes(fill = IL1B), alpha = 0.4, size = 3, colour="transparent") +
#   stat_ellipse(
#     aes(
#       # fill = paste(Infection, Status)
#       colour = Infection,
#       linetype = Status
#     ),
#     subset(il1b.data, Time == "6h"),
#     alpha = 0.4,
#     level = 0.5,
#     # geom = "polygon",
#     geom = "path"
#   ) +
#   scale_fill_gradient() +
#   theme_bw()

ggplot(il1b.data, aes(tSNE1, tSNE2, shape = paste(Infection, Status))) +
  geom_rug() +
  stat_ellipse(
    aes(
      fill = paste(Infection, Status)
    ),
    subset(il1b.data, Time == "6h"),
    alpha = 0.4,
    level = 0.5,
    geom = "polygon"
    # geom = "path"
  ) +
  geom_point(aes(colour = IL1B), alpha = 0.75, size = 1.5) +
  scale_color_gradient2(
    low = col.gradient[11], mid = col.gradient[6], high = col.gradient[1],
    midpoint = max(norm_exprs(sce.endo)) / 2 - min(norm_exprs(sce.endo)) / 2
  ) +
  scale_fill_manual(values = col.treatment) +
  labs(fill = "Treatment", shape = "Treatment") +
  theme_bw()
ggsave("28_out/tSNE_ellipses_themeBW.pdf", height = 5, width = 8)

ggplot(il1b.data, aes(tSNE1, tSNE2, shape = paste(Infection, Status))) +
  # geom_rug() +
  # stat_ellipse(
  #   aes(
  #     fill = paste(Infection, Status)
  #   ),
  #   subset(il1b.data, Time == "6h"),
  #   alpha = 0.4,
  #   level = 0.5,
  #   geom = "polygon"
  #   # geom = "path"
  # ) +
  geom_point(aes(colour = IL1B), alpha = 0.75, size = 1.5) +
  scale_color_gradient2(
    low = col.gradient[11], mid = col.gradient[6], high = col.gradient[1],
    midpoint = max(norm_exprs(sce.endo)) / 2 - min(norm_exprs(sce.endo)) / 2
  ) +
  scale_fill_manual(values = col.treatment) +
  labs(fill = "Treatment", shape = "Treatment") +
  theme_bw()
ggsave("28_out/tSNE_themeBW.pdf", height = 5, width = 8)

ggplot(il1b.data, aes(tSNE1, tSNE2, shape = paste(Infection, Status))) +
  # geom_rug() +
  # stat_ellipse(
  #   aes(
  #     fill = paste(Infection, Status)
  #   ),
  #   subset(il1b.data, Time == "6h"),
  #   alpha = 0.4,
  #   level = 0.5,
  #   geom = "polygon"
  #   # geom = "path"
  # ) +
  geom_point(aes(colour = IL1B), alpha = 0.75, size = 2) +
  scale_color_gradient2(
    low = col.gradient[11], mid = col.gradient[6], high = col.gradient[1],
    midpoint = max(norm_exprs(sce.endo)) / 2 - min(norm_exprs(sce.endo)) / 2
  ) +
  scale_fill_manual(values = col.treatment) +
  labs(fill = "Treatment", shape = "Treatment") +
  theme_grey()
ggsave("28_out/tSNE_themeGrey.pdf", height = 5, width = 8)

ggplot(il1b.data, aes(tSNE1, tSNE2, shape = paste(Infection, Status))) +
  # geom_rug() +
  # stat_ellipse(
  #   aes(
  #     fill = paste(Infection, Status)
  #   ),
  #   subset(il1b.data, Time == "6h"),
  #   alpha = 0.4,
  #   level = 0.5,
  #   geom = "polygon"
  #   # geom = "path"
  # ) +
geom_point(aes(colour = IL1B), alpha = 0.75, size = 2) +
  scale_color_gradient2(
    low = col.gradient[11], mid = col.gradient[6], high = col.gradient[1],
    midpoint = max(norm_exprs(sce.endo)) / 2 - min(norm_exprs(sce.endo)) / 2
  ) +
  scale_fill_manual(values = col.treatment) +
  labs(fill = "Treatment", shape = "Treatment") +
  theme_dark()
ggsave("28_out/tSNE_themeDark.pdf", height = 5, width = 8)

ggplot(il1b.data, aes(tSNE1, tSNE2, shape = paste(Infection, Status))) +
  # geom_rug() +
  # stat_ellipse(
  #   aes(
  #     fill = paste(Infection, Status)
  #   ),
  #   subset(il1b.data, Time == "6h"),
  #   alpha = 0.4,
  #   level = 0.5,
  #   geom = "polygon"
  #   # geom = "path"
  # ) +
geom_point(aes(colour = IL1B), alpha = 0.75, size = 1.5) +
  scale_color_gradient2(
    low = col.gradient[11], mid = col.gradient[5], high = col.gradient[1],
    midpoint = max(norm_exprs(sce.endo)) / 2 - min(norm_exprs(sce.endo)) / 2
  ) +
  scale_fill_manual(values = col.treatment) +
  labs(fill = "Treatment", shape = "Treatment") +
  theme_bw()
ggsave("28_out/tSNE_themeBW_midOrange.pdf", height = 5, width = 8)
