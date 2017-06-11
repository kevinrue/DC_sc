library(readxl)
library(reshape2)

qPCR <- read_excel("~/git/DC_sc/06_qPCR/qPCR.xlsx", sheet = "Sheet1")

qPCR[qPCR == 0] <- NA

repNames <- paste("Replicate", 1:4)

qPCR$Gene <- factor(qPCR$Gene, unique(qPCR$Gene))
qPCR$Status <- factor(qPCR$Status, unique(qPCR$Status))

colours.infection <- brewer.pal(3, "Set1")[2:3]
names(colours.infection) <- c("STM-D23580", "STM-LT2")

qPCR$logMean <- rowMeans(log10(qPCR[,repNames]), na.rm = TRUE)
qPCR$logMeanSEM <- apply(log10(qPCR[,repNames]), 1, function(x){
  sd(x, na.rm = TRUE)/sqrt(length(x))})
qPCR$logMeanMinusSEM <- qPCR$logMean - qPCR$logMeanSEM
qPCR$logMeanPlusSEM <- qPCR$logMean + qPCR$logMeanSEM
dim(qPCR)

ggplot(qPCR, aes(Infection, 10^logMean)) +
  facet_grid(Gene ~ Status, scales = "free_y") +
  geom_col(aes(fill = Infection), position = "dodge") + # 10^
  geom_errorbar(aes(ymin = 10^logMeanMinusSEM, ymax = 10^logMeanPlusSEM), width = 0.2) +
  scale_fill_manual(values = colours.infection) +
  # scale_y_log10() +
  theme_minimal()

res.test <- data.frame()
for (x_gene in levels(qPCR$Gene)){
  for (x_status in levels(qPCR$Status)){
    message(sprintf("%s %s", x_gene, x_status))
    d23 <- log10(as.numeric(subset(
      qPCR,
      Gene == x_gene & Status == x_status & Infection == "STM-D23580",
      repNames)))
    lt2 <- log10(as.numeric(subset(
      qPCR,
      Gene == x_gene & Status == x_status & Infection == "STM-LT2",
      repNames)))
    t.tmp <- tidy(t.test(d23, lt2))
    label.y <- 10^max(subset(
      qPCR, Gene == x_gene & Status == x_status, "logMeanPlusSEM"), na.rm = TRUE)
    t.tmp <- cbind(
      Gene = x_gene,
      Status = x_status,
      LabelY = label.y,
      t.tmp
    )
    res.test <- rbind(res.test, t.tmp)
  }
}
View(res.test)

res.test$sig <- as.character(symnum(
  res.test$p.value, corr = FALSE,
  cutpoints = c(0,  .001,.01,.05, .1, 1),
  symbols = c("***","**","*",""," "), na = " "))

gg <- ggplot(as.data.frame(qPCR), aes(Status, 10^logMean)) +
  facet_grid(Gene ~ ., scales = "free_y") +
  geom_col(aes(fill = Infection, colour = Infection), position = "dodge") + # 10^
  geom_errorbar(
    aes(ymin = 10^logMeanMinusSEM, ymax = 10^logMeanPlusSEM, colour = Infection),
    width = 0.15, size = 0.5, position = position_dodge(0.9)
  ) +
  geom_text(aes(y = LabelY, label = sig), res.test) +
  scale_fill_manual(values = colours.infection) +
  scale_colour_manual(values = colours.infection) +
  labs(y = "Fold-change") + theme_minimal() +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

gg

for (x_gene in levels(qPCR$Gene)){
  for (x_status in levels(qPCR$Status)){
    x_name <- sprintf("%s - %s", x_gene, x_status)
    f_name <- sprintf("%s_%s.pdf", x_gene, x_status)
    message(x_name)
    ggData <- subset(qPCR, Gene == x_gene & Status == x_status)
    ggText <- subset(res.test, Gene == x_gene & Status == x_status)
    gg <- ggplot(as.data.frame(ggData), aes(Infection, 10^logMean)) +
      geom_col(aes(fill = Infection, colour = Infection)) +
      geom_errorbar(
        aes(ymin = 10^logMeanMinusSEM, ymax = 10^logMeanPlusSEM, colour = Infection),
        width = 0.15, size = 0.5, position = position_dodge(0.9)
      ) +
      scale_fill_manual(values = colours.infection) +
      scale_colour_manual(values = colours.infection) +
      labs(y = "Fold-change") + theme_minimal() +
      theme(
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold")
      ) + labs(title = x_name)
    ggsave(f_name, gg, height = 4, width = 6)
  }
}
