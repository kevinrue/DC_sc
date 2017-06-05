library(readxl)
library(broom)
library(RColorBrewer)
library(ggplot2)

uptake <- read_excel(
  "~/git/DC_sc/05_bacterial_growth/DC_bacterial_growth.xlsx",
  sheet = "Uptake"
)

repNames <- paste("Replicate", 1:3)

uptake$Mean <- rowMeans(uptake[,repNames])
uptake$MeanSEM <- apply(uptake[,repNames], 1, function(x){sd(x)/sqrt(length(x))})
uptake$MeanMinusSEM <- uptake$Mean - uptake$MeanSEM
uptake$MeanPlusSEM <- uptake$Mean + uptake$MeanSEM
dim(uptake)

colours.infection <- brewer.pal(3, "Set1")[2:3]
names(colours.infection) <- c("STM-D23580", "STM-LT2")

ggplot(uptake, aes(Infection)) +
  geom_point(aes(y=Mean), colour = "black") +
  geom_point(aes(y=MeanMinusSEM), colour = "green") +
  geom_point(aes(y=MeanPlusSEM), colour = "red")

uptake.test <- data.frame()
d23 <- as.numeric(subset(
  uptake,
  Infection == "STM-D23580",
  repNames))
lt2 <- as.numeric(subset(
  uptake,
  Infection == "STM-LT2",
  repNames))
t.tmp <- tidy(t.test(d23, lt2))
t.tmp <- cbind(
  Max = max(c(d23, lt2)),
  t.tmp
)
uptake.test <- rbind(uptake.test, t.tmp)

uptake.test$sig <- symnum(
  uptake.test$p.value, corr = FALSE,
  cutpoints = c(0,  .001,.01,.05, .1, 1),
  symbols = c("***","**","*",""," "), na = " ")

gg <- ggplot(
  as.data.frame(uptake),
  aes(Infection, Mean, fill = Infection, colour = Infection)
) +
  geom_errorbar(
    aes(ymin = MeanMinusSEM, ymax = MeanPlusSEM),
    width = 0.15, size = 0.5
  ) +
  geom_col(aes(group = Infection)) +
  scale_fill_manual(values = colours.infection) +
  scale_colour_manual(values = colours.infection) +
  theme_minimal() + labs(y = "% uptake") +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

ggsave("uptake_bars.pdf", gg, height = 5, width = 6)
