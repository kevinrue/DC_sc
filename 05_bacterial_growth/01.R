library(readxl)
library(broom)
library(RColorBrewer)
library(ggplot2)

cfuml <- read_excel("~/git/DC_sc/05_bacterial_growth/DC_bacterial_growth.xlsx")

repNames <- paste("Replicate", 1:3)

cfuml$Mean <- rowMeans(cfuml[,repNames])
cfuml$MeanSEM <- apply(cfuml[,repNames], 1, function(x){sd(x)/sqrt(length(x))})
cfuml$MeanMinusSEM <- cfuml$Mean - cfuml$MeanSEM
cfuml$MeanPlusSEM <- cfuml$Mean + cfuml$MeanSEM
cfuml$logMean <- rowMeans(log10(cfuml[,repNames]))
cfuml$logMeanSEM <- apply(log10(cfuml[,repNames]), 1, function(x){sd(x)/sqrt(length(x))})
dim(cfuml)

cfuml$Time <- factor(cfuml$Time, c("45min", "75min", "2h", "4h", "6h"))

colours.infection <- brewer.pal(3, "Set1")[2:3]
names(colours.infection) <- c("STM-D23580", "STM-LT2")

ggplot(cfuml, aes(Time)) + facet_grid(Infection ~ .) +
  geom_point(aes(y=Mean), colour = "black") +
  geom_point(aes(y=MeanMinusSEM), colour = "green") +
  geom_point(aes(y=MeanPlusSEM), colour = "red") + scale_y_log10()

cfu.test <- data.frame()
for (x_time in levels(cfuml$Time)){
  message(sprintf("%s", x_time))
  d23 <- as.numeric(subset(
    cfuml,
    Time == x_time & Infection == "STM-D23580",
    repNames))
  lt2 <- as.numeric(subset(
    cfuml,
    Time == x_time & Infection == "STM-LT2",
    repNames))
  t.tmp <- tidy(t.test(d23, lt2))
  t.tmp <- cbind(
    Time = x_time,
    Max = max(c(d23, lt2)),
    t.tmp
  )
  cfu.test <- rbind(cfu.test, t.tmp)
}

cfu.test$sig <- symnum(
  cfu.test$p.value, corr = FALSE,
  cutpoints = c(0,  .001,.01,.05, .1, 1),
  symbols = c("***","**","*",""," "), na = " ")

gg <- ggplot(
  as.data.frame(cfuml),
  aes(Time, Mean, fill = Infection, colour = Infection)
) +
  geom_errorbar(
    aes(ymin = MeanMinusSEM, ymax = MeanPlusSEM),
    width = 0.15, size = 0.5
  ) +
  geom_line(aes(group = Infection)) +
  scale_fill_manual(values = colours.infection) +
  scale_y_log10(
    breaks = 10^(4:8),
    limits = c(10^4, 10^8)
  ) +  annotation_logticks(sides = "l") +
  scale_colour_manual(values = colours.infection) +
  labs(y = "log10 (CFU / mL)") + theme_minimal() +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

ggsave("CFUml_bars.pdf", gg, height = 5, width = 6)
