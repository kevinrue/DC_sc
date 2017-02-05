
groupName = "4h_LT2_exposed"
cellIndex <- which(sce.norm$Group == groupName)
length(cellIndex)

CountsQC = counts(sce.norm)[,cellIndex]
dim(CountsQC)
typeof(CountsQC)
class(CountsQC)

TechQC = grepl("^ERCC-", rownames(CountsQC))
table(TechQC)
sum(tail(TechQC, n=39)) # ERCC are in the last rows of the count matrix

SpikeInfo <- mixC[, c("ID","ERCCmolecules")]
colnames(SpikeInfo) <- c("Name", "molecules_in_each_chamber")
head(SpikeInfo, n=4)

SpikeInfoQC <- SpikeInfo[SpikeInfo$Name %in% rownames(CountsQC)[TechQC],]

Data = newBASiCS_Data(Counts = CountsQC,
                      Tech = TechQC,
                      SpikeInfo = SpikeInfoQC,
                      BatchInfo = sce.norm$Plate[cellIndex])

MCMC_Output <- BASiCS_MCMC(Data,
                           N = 20,
                           Thin = 2,
                           Burn = 4)

N = 200; Thin = 10; Burn = 100 # published N = 20000; Thin = 10; Burn = 10000
MCMC_Output <- BASiCS_MCMC(
  Data, N = N, Thin = Thin, Burn = Burn,
  PriorDelta = 'log-normal')
