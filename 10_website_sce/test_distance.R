sce.sc <- readRDS("rds/sce.sc.SCE.rds")

assayNames(sce.sc)
class(assays(sce.sc))

class(assay(sce.sc, "counts"))

counts <- assay(sce.sc, "counts")
dim(counts)

numDetect <- rowSums(counts > 0)

hist(numDetect)

counts <- counts[numDetect > 0,]

distanceToCell <- function(matrix, cell){
  refCell <- matrix[,cell]
  otherCells <- matrix[,-cell]
  stopifnot(ncol(otherCells) == (ncol(matrix)-1))
  otherCells <- abs(otherCells - refCell)
  return(colSums(otherCells))
}

colData(sce.sc)[1, c("Time","Infection","Status")]

d1 <- distanceToCell(counts, 1)
plot(sort(d1))

d1 <- data.frame(
  Distance = distanceToCell(counts, 1),
  colData(sce.sc)[-1, c("Time","Infection","Status")]
)
d1$Rank <- rank(d1$Distance)

ggplot(d1) +
  geom_point(aes(Rank, log10(Distance), colour=Time, shape=Status), size=4, alpha=1/4)

logDistanceToCell <- function(matrix, cell){
  refCell <- matrix[,cell]
  otherCells <- matrix[,-cell]
  stopifnot(ncol(otherCells) == (ncol(matrix)-1))
  otherCells <- log10(1 + abs(otherCells - refCell))
  return(colSums(otherCells))
}

dlog1 <- data.frame(
  log10Distance = logDistanceToCell(counts, 1),
  colData(sce.sc)[-1, c("Time","Infection","Status")]
)
dlog1$Rank <- rank(dlog1$log10Distance)

ggplot(dlog1) +
  geom_point(aes(Rank, log10Distance, colour=Time, shape=Status)) +
  labs(
    subtitle = with(
      data.frame(colData(sce.sc)[1,]),
      sprintf("%s : %s | %s | %s", Sample, Time, Infection, Status)
    )
  )

plotManhattanDistanceToCell <- function(counts, index, Ylabel="counts"){
  distanceDataFrame <- data.frame(
    Distance = distanceToCell(counts, index),
    colData(sce.sc)[-index, c("Time","Infection","Status")] # hard-code
  )
  distanceDataFrame$Rank <- rank(distanceDataFrame$Distance)
  ggplot(distanceDataFrame) +
    geom_point(aes(Rank, Distance / 1E6, colour=Time, shape=Status)) +
    labs(
      title = expression(log[10]*" Distance"),
      subtitle = with(
        data.frame(colData(sce.sc)[index,]),
        sprintf("[%i] %s : %s | %s | %s", index, Sample, Time, Infection, Status)
      ),
      y = sprintf("Manhattan %s distance (million)", Ylabel),
      x = "Distance rank"
    )
}

plotManhattanDistanceToCell(counts, sample(ncol(counts), 1))

# Note: depth difference between cells is not considered so far
# e.g. use CPM to normalise between cells

cellSize <- colSums(counts)
cpm <- t(t(counts) / cellSize)

plotManhattanDistanceToCell(cpm, sample(ncol(cpm), 1), "CPM")
