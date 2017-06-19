library(readxl)
cytokines <- read_excel("cytokines complete.xlsx")

repNames <- paste("Replicate", 1:3)

cytokines$Cytokine <- factor(cytokines$Cytokine, c("IL12","IL1b","TNFa","IL6","IL10"))
cytokines$Mean <- rowMeans(cytokines[,repNames])
cytokines$SEM <- apply(cytokines[,repNames], 1, function(x){sd(x)/sqrt(length(x))})
dim(cytokines)

cyto.1 <- subset(cytokines, Cytokine == "IL10")
dim(cyto.1)

cyto.2 <- subset(cytokines, Cytokine != "IL10")
dim(cyto.2)

library(broom)
t.df <- data.frame()
for (x_cytokine in unique(cytokines$Cytokine)){
  for (x_moi in sort(unique(cytokines$MOI))){
    for (x_time in unique(cytokines$Time)){
      message(sprintf("%s %s %s", x_cytokine, x_moi, x_time))
      d23 <- as.numeric(subset(
        cytokines,
        Cytokine == x_cytokine & MOI == x_moi & Time == x_time & Infection == "STM-D23580",
        repNames))
      lt2 <- as.numeric(subset(
        cytokines,
        Cytokine == x_cytokine & MOI == x_moi & Time == x_time & Infection == "STM-LT2",
        repNames))
      t.tmp <- tidy(t.test(d23, lt2))
      t.tmp <- cbind(
        Cytokine = x_cytokine, MOI = x_moi, Time = x_time,
        Max = max(c(d23, lt2)),
        t.tmp
      )
      t.df <- rbind(t.df, t.tmp)
    }
  }
}

t.df$sig <- symnum(
  t.df$p.value, corr = FALSE,
  cutpoints = c(0,  .001,.01,.05, .1, 1),
  symbols = c("***","**","*",""," "), na = " ")

t.1 <- subset(t.df, Cytokine == "IL10")
dim(t.1)

t.2 <- subset(t.df, Cytokine != "IL10")
dim(t.2)

library(RColorBrewer)
colours.infection <- brewer.pal(3, "Set1")[2:3]
names(colours.infection) <- rev(unique(moi$Infection))

library(ggplot2)
basePlot.1 <-
  ggplot(cyto.1, aes(Time, Mean)) +
  geom_errorbar(
    aes(ymin = Mean - SEM, ymax = Mean + SEM, colour = Infection),
    position = "dodge", width = 0.9, size = 0.5
  ) +
  geom_col(
    aes(fill = Infection, colour = Infection),
    position = "dodge") +
  geom_text(aes(Time, Max*1.1, label = sig), t.1) +
  facet_grid(Cytokine ~ MOI, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = colours.infection) +
  scale_colour_manual(values = colours.infection) +
  labs(y = "pg/mL") + theme_bw()

ggsave("TheUglyIL10Cytokine.pdf", basePlot.1, width = 8, height = 2)


basePlot.2 <-
  ggplot(cyto.2, aes(Time, Mean)) +
  geom_errorbar(
    aes(ymin = Mean - SEM, ymax = Mean + SEM, colour = Infection),
    position = "dodge", width = 0.9, size = 0.5
  ) +
  geom_col(
    aes(fill = Infection, colour = Infection),
    position = "dodge") +
  geom_text(aes(Time, Max*1.1, label = sig), t.2) +
  facet_grid(Cytokine ~ MOI, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = colours.infection) +
  scale_colour_manual(values = colours.infection) +
  labs(y = "pg/mL") + theme_bw() + theme(legend.position = "bottom")

ggsave("TheOtherCuteCytokines.pdf", basePlot.2, width = 8, height = 6)
