library(readxl)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

cfu <- read_excel(
  "~/git/DC_sc/05_bacterial_growth/cfu stacked plot.xlsx", sheet = "Clean"
)
View(cfu)

cfuLevels <- c("0" ,"1-3", "4-7", "> 8")
cfu$CFU <- factor(as.character(cfu$CFU), rev(cfuLevels))

cfu$Percentage <- apply(cfu, 1, function(x){
  status = x[2]
  infection = x[3]
  time = x[4]
  count = as.numeric(x[5])
  total = sum(subset(
    cfu,
    Status == status & Infection == infection & Time == time,
    "Count"))
  perc = count / total
  perc
})

display.brewer.all()

# gradient.cfu <- brewer.pal(9, "Reds")[c(2,4,6,8)]
# gradient.cfu <- brewer.pal(9, "Blues")[c(2,4,6,8)]
# gradient.cfu <- brewer.pal(9, "BuPu")[c(1,3,5,7)]
gradient.cfu <- brewer.pal(9, "RdGy")[c(7,4,3,2)]
# gradient.cfu <- c(
#   brewer.pal(9, "RdGy")[7],
#   brewer.pal(9, "Reds")[c(2,4,6)]
# )

names(gradient.cfu) <- cfuLevels

brks <- c(0, 0.25, 0.5, 0.75, 1)

gg <- ggplot(cfu) +
  geom_col(aes(Time, Percentage, fill = CFU)) +
  facet_grid(~ Status + Infection) +
  scale_fill_manual(values = rev(gradient.cfu)) +
  scale_y_continuous(breaks = brks, labels = scales::percent(brks)) +
  theme_minimal()
gg

ggsave("Percentage_CFUlevel.pdf", gg, width = 6, height = 4)

