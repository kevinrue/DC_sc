sce.norm <- readRDS("rds/sce.norm.tSNE.rds")

varLabels(sce.norm)

tmpGG <- data.frame(
  reducedDimension(sce.norm),
  Features = sce.norm$total_features,
  Time = sce.norm$Time,
  Treatment = gsub("_", " ", sce.norm$Treatment)
)

ggplot(tmpGG, aes(X1, X2, shape = Time, colour = Features)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(x="Dimension 1",y="Dimension 2",title="2h")

ggsave("17_out/tSNE_total_features.pdf", height = 4, width = 5)
