if(!exists('doneSetup'))
  source('00-setup.R')

# FEL2008 ----
D = read.table('../data/FEL2008_Data.csv',
               sep = ",", header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)
systems = D[,1]

## Full dataset ----
C  = D[,2]
uC = D[,3]
R  = D[,4]
uR = D[,5] / 1.96
E  = R-C
uE = sqrt(uC^2+uR^2)
N = length(E)

## ELiminate outlier (SiH) ----
io = which(grepl('SiH',systems))
D0 = D[-io,]

S0  = D0[,1]
C0  = D0[,2]
uC0 = D0[,3]  # Hyp: std. unc.
R0  = D0[,4]
uR0 = D0[,5] / 1.96
E0  = R0-C0
uE0 = sqrt(uC0^2+uR0^2)
z   = E0/uE0
N0  = length(z)


label = 1
png(file = paste0(figDir,'/dist_FEL2008.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  C, E, uy = uE,
  plotGauss = TRUE, nclass=17,
  xlab = 'Calc. atom. energ. [kcal/mol]',
  ylab = 'Error [kcal/mol]',
  plotBA = TRUE,
  outLiers = TRUE, p=0.95,
  labels = systems,
  xlim = c(-150,250),
  ylim = c(-15,20),
  topMar = 3,
  gPars = gPars,
  scalePoints = 0.25
)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()

label = 2
png(file = paste0(figDir,'/dist_FEL2008_zscores.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  uE0, z,
  logX = TRUE,
  plotGauss = TRUE, nclass = 17,
  xlab = 'Pred. unc. [kcal/mol]',
  ylab = 'z-scores',
  plotBA = TRUE,
  outLiers = TRUE, p=0.95,
  # xlim = c(-150,250),
  # ylim = c(-1.25,1.25),
  topMar = 3,
  labels = S0,
  gPars = gPars,
  scalePoints = 0.25
)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()

label =3
png(file = paste0(figDir,'/LZV_FEL2008.png'),
    width = gPars$reso, height = gPars$reso)
res = ErrViewLib::plotLZV(
  uE0, z,
  logX = TRUE,
  xlab = 'Pred. unc. [kcal/mol]',
  # xlim = c(0.02,1),
  # ylim      = c(0,1.2),
  label     = label,
  gPars     = gPars
)
dev.off()
