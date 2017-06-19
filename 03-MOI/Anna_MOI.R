library(readxl)
MOI <- read_excel("MOI_complete.xlsx")
View(MOI)

library(RColorBrewer)
annasFavoriteColours <- brewer.pal(3, "Set1")[2:3]
names(annasFavoriteColours) <- rev(unique(MOI$Infection))

library(ggplot2)
basePlot <- ggplot(MOI) +
  geom_col(
    aes(Time, Value, fill = Infection),
    position = "dodge") +
  facet_grid(~ MOI) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = annasFavoriteColours) +
  labs(y = "CFU/mL") +
  scale_y_continuous(
    breaks = seq(0, max(MOI$Value), 10E3),
    labels = c(0, sprintf("%i,000", seq(10E3, max(MOI$Value), 10E3)/1E3)))

cutePlotWithAnnasColours <- basePlot + theme_bw()
# basePlot + theme_minimal() + labs(title = "minimal")

ggsave("MyLittleMOIplot.pdf", cutePlotWithAnnasColours, width = 8, height = 2)
