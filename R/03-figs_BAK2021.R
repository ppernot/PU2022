if(!exists('doneSetup'))
  source('00-setup.R')

# Stoechiometry functions ----
elements    = c("H", "C", "N", "O", "F")
filterFormula <- function(sp) {
  # Filter out non-stoechimetric notations
  # to enable composition calculation in get.atoms()
  sp1 = sub("\\..*$", "", sp)
  return(sp1)
}
calcAtoms <- function(formula) {
  # Transform filtered formula to composition vector
  atoms = CHNOSZ::i2A(formula)
  compo = matrix(0, nrow = 1, ncol = length(elements))
  colnames(compo) = elements
  compo[1, colnames(atoms)] = atoms
  return(compo[1, ])
}
get.atoms <- function(sp) {
  # Transform formula to composition vector
  sp1 = filterFormula(sp)
  tryCatch(calcAtoms(sp1), error = function(x) rep(NA, length(elements)))
}


# BAK2021 ----
D = read.table('../data/BAK2021_Data.csv',
               sep = ",", header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)
systems = D[,1]

compo = matrix(0, nrow = length(systems), ncol = length(elements))
colnames(compo) = elements
heteroFrac = c()
for (i in seq_along(systems)) {
  atoms = get.atoms(systems[i])
  compo[i, names(atoms)] = atoms
  heteroFrac[i] = (compo[i,'N']+compo[i,'O']+compo[i,'F'])/sum(compo[i,])
}

## Full dataset ----

R  = D[,2]
C  = D[,3]
uC = D[,4]/1.96 # Convert u95 to u for consistency in plot functions
E  = R-C
N = length(E)

# 95 % PICP
# Bootstrap uncertainty on I95 coverage
df = data.frame(E=E, Up=D[,4])
bootPICP(df, numDig = 1)


label = 1
png(file = paste0(figDir,'/dist_BAK2021.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  C, E, uy = uC,
  plotGauss = TRUE, nclass=17,
  xlab = 'Calc. ZPE [kcal/mol]',
  ylab = 'Error [kcal/mol]',
  ylim = c(-1,1),
  labels = labels,
  topMar = 3,
  scalePoints = 0.25,
  gPars = gPars
)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()

label = label + 1
png(file = paste0(figDir,'/pcoverage_BAK2021.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotPcoverage(
  as.data.frame(cbind(R = R, D = C, uP = uC)),
  mycols    = c(2,5,7),
  prob      = 0.95,
  nBin      = 2,
  slide     = TRUE,
  xlab      = 'Calc. ZPE [kcal/mol]',
  # xlab      = 'Pred. unc. [kcal/mol]',
  ylim      = c(0.7,1),
  label = label,
  gPars     = gPars
)
dev.off()

label = label + 1
png(file = paste0(figDir,'/pcoverage_UP_BAK2021.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotPcoverage(
  as.data.frame(cbind(R = R, D = C, uP = uC)),
  ordX      = uC,
  mycols    = c(2,5,7),
  nBin = 2, slide = TRUE,
  prob      = 0.95,
  xlab      = 'Pred. unc. [kcal/mol]',
  ylim      = c(0.7,1),
  label = label,
  gPars     = gPars
)
dev.off()

label = label + 1
png(file = paste0(figDir,'/pcoverage_Het_BAK2021.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotPcoverage(
  as.data.frame(cbind(R = R, D = C, uP = uC)),
  ordX      = heteroFrac,
  nBin      = 2,
  slide     = TRUE,
  mycols    = c(2,5,7),
  prob      = 0.95,
  xlab      = 'Fraction of hereroatoms',
  ylim      = c(0.7,1),
  label = label,
  gPars     = gPars
)
dev.off()
