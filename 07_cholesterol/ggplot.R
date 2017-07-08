library(readxl)
cholesterol <- read_excel("~/git/DC_sc/07_cholesterol/clean_data.xlsx")
colnames(cholesterol)

cholesterol$Compound <- factor(cholesterol$Compound, unique(cholesterol$Compound))

cholesterol$Infection <- factor(cholesterol$Infection, unique(cholesterol$Infection))

library(ggplot2)

ggplot(
  subset(cholesterol, Compound != "Ester"),
  aes(
    factor(gsub("-", "\n", Infection), gsub("-", "\n", levels(Infection))),
    Value)
  ) +
  facet_wrap(~ Compound, ncol = 1, scales = "free_y") +
  geom_violin(
    aes(fill = Compound),
    draw_quantiles = c(0.25, 0.5, 0.75),
    linetype = "dashed", colour = "grey"
  ) +
  labs(x = "Infection") +
  theme_bw()
  # geom_violin(aes(Infection, Value), draw_quantiles = c(0.25, 0.5, 0.75))

ggsave("violin_(noEster).pdf", height = 8, width = 6)
