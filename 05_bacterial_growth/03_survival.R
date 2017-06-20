library(readxl)
library(broom)
library(RColorBrewer)
library(ggplot2)

survival <- read_excel("DC_bacterial_growth.xlsx",sheet = "Survival")

repNames <- paste("Replicate", 1:3)

survival$Time <- factor(survival$Time, c("75min", "2h", "4h", "6h"))

survival$Mean <- rowMeans(survival[,repNames])
survival$MeanSEM <- apply(survival[,repNames], 1, function(x){sd(x)/sqrt(length(x))})
survival$MeanMinusSEM <- survival$Mean - survival$MeanSEM
survival$MeanPlusSEM <- survival$Mean + survival$MeanSEM
dim(survival)

colours.infection <- brewer.pal(3, "Set1")[2:3]
names(colours.infection) <- c("STM-D23580", "STM-LT2")

ggplot(survival, aes(Time)) + facet_grid(~ Infection) +
  geom_point(aes(y=Mean), colour = "black") +
  geom_point(aes(y=MeanMinusSEM), colour = "green") +
  geom_point(aes(y=MeanPlusSEM), colour = "red")

survival.test <- data.frame()
for (x_time in levels(survival$Time)){
  message(sprintf("%s", x_time))
  d23 <- as.numeric(subset(
    survival,
    Time == x_time & Infection == "STM-D23580",
    repNames))
  lt2 <- as.numeric(subset(
    survival,
    Time == x_time & Infection == "STM-LT2",
    repNames))
  t.tmp <- tidy(t.test(d23, lt2))
  t.tmp <- cbind(
    Time = x_time,
    Max = max(c(d23, lt2)),
    t.tmp
  )
  survival.test <- rbind(survival.test, t.tmp)
}

survival.test$sig <- symnum(
  survival.test$p.value, corr = FALSE,
  cutpoints = c(0,  .001,.01,.05, .1, 1),
  symbols = c("***","**","*",""," "), na = " ")

gg <- ggplot(
  as.data.frame(survival),
  aes(Time, Mean, fill = Infection, colour = Infection)
  ) +
  geom_errorbar(
    aes(ymin = MeanMinusSEM, ymax = MeanPlusSEM),
    position = position_dodge(width = 0.9), width = 0.25, size = 0.5
  ) +
  geom_col(aes(group = Infection), position = "dodge") +
  scale_fill_manual(values = colours.infection) +
  scale_colour_manual(values = colours.infection) +
  theme_minimal() + labs(y = "% survival") +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

ggsave("survival_bars.pdf", gg, height = 5, width = 6)

gg <- ggplot(
  as.data.frame(survival),
  aes(Time, Mean, fill = Infection, colour = Infection)
) +
  geom_errorbar(
    aes(ymin = MeanMinusSEM, ymax = MeanPlusSEM),
    position = position_dodge(width = 0.9), width = 0.25, size = 0.5
  ) +
  geom_col(aes(group = Infection), position = "dodge") +
  scale_fill_manual(values = colours.infection) +
  scale_colour_manual(values = colours.infection) +
  theme_minimal() + labs(y = "% survival") +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("survival_bars_legendBottom.pdf", gg, height = 5, width = 6)
