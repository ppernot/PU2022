if(!exists('doneSetup'))
  source('00-setup.R')

gParsPlot = ErrViewLib::setgPars()

# PRO2021 ----
D = read.table('../data/PRO2021_Data.csv',
               sep = ",", header = TRUE,
               check.names = FALSE,
               stringsAsFactors = FALSE)

R  = D[,1]
C  = D[,2]
ua = D[,3] # U95
ub = D[,4] # U95
E  = R-C
N = length(E)

df = data.frame(E=E,Up=D[,3])
print(bootPICP(df, numDig = 1))

df = data.frame(E=E,Up=D[,4])
print(bootPICP(df, numDig = 1))

label = 1
png(file = paste0(figDir,'/dist_PRO2021a.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  C, E, ua/1.96, # Compensate for expansion to 2-sigma in plotDistHist
  plotGauss = TRUE, nclass=17,
  xlab = 'Log calculated rate',
  ylab = 'Error',
  plotBA = TRUE,
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
png(file = paste0(figDir,'/LZV_PRO2021a.png'),
    width = gPars$reso, height = gPars$reso)
res = ErrViewLib::plotLZV(
  ua/1.96, E/(ua/1.96),
  logX = TRUE,
  xlab      = 'Prediction uncertainty',
  xlim = c(0.02,1),
  ylim      = c(0,1.2),
  label     = label,
  gPars     = gPars
)
dev.off()

label = 3
png(file = paste0(figDir,'/LZV_PRO2021b.png'),
    width = gPars$reso, height = gPars$reso)
res = ErrViewLib::plotLZV(
  ub/1.96, E/(ub/1.96),
  logX = TRUE,
  xlab      = 'Prediction uncertainty',
  xlim = c(0.02,1),
  ylim      = c(0,1.2),
  label     = label,
  gPars     = gPars
)
dev.off()

# For reference ####

png(file = paste0(figDir,'/dist_PRO2021b.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  C, E, ub/1.96,
  plotGauss = TRUE,
  nclass=17,
  xlab = 'Log calculated rate',
  ylab = 'Error',
  plotBA = TRUE,
  topMar = 3,
  gPars = gPars,
  scalePoints = 0.25
)
dev.off()

png(file = paste0(figDir,'/LCP_PRO2021a.png'),
    width = gPars$reso, height = gPars$reso)
res = ErrViewLib::plotPcoverage(
  as.data.frame(cbind(D=C, R=R, uP=ua/1.96)),
  mycols    = c(2,5,7),
  prob      = 0.95,
  xlab      = 'Log caculated rate',
  ylim      = c(0.94,1.0),
  gPars     = gPars
)
dev.off()

png(file = paste0(figDir,'/LCP_PRO2021b.png'),
    width = gPars$reso, height = gPars$reso)
res = ErrViewLib::plotPcoverage(
  as.data.frame(cbind(D=C, R=R, uP=ub/1.96)),
  mycols    = c(2,5,7),
  prob      = 0.95,
  xlab      = 'Log caculated rate',
  ylim      = c(0.94,1.0),
  gPars     = gPars
)
dev.off()

