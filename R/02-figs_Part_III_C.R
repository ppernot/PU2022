if(!exists('doneSetup'))
  source('00-setup.R')

set.seed(123)

N     = 1000
V     = rchisq(N, df=4)
uE    = 0.01 * sqrt(V / mean(V))
uE0   = uE
calc  = runif(N)
nu    = 5
err   = uE * rt(N, df=nu) * sqrt(nu-2)/sqrt(nu)
ref   = calc + err
z     = err / uE

# Stats
print(ErrViewLib::kurtcs(err))
print(mean(err))
print(sd(err))
print(ErrViewLib::kurtcs(z))



# Calibrated set ####

# PICP
picp = mean( z >= -1.96 & z <= 1.96 )
print(picp)
ci = DescTools::BinomCI(picp * N, N, method = BPCImethod)[,2:3]
print(ci)
sha = sqrt(mean(uE^2))
print(sha)


subset = 1:250
label = 1
png(file = paste0(figDir,'/distErrorsNew.png'),
    width= gPars$reso,height=gPars$reso)
ErrViewLib::plotDistHist(
  calc[subset], err[subset], uE[subset],
  main = '',
  plotGauss = TRUE,
  plotReg = FALSE, degree = 1,
  xlab = 'Calculated value [a.u.]',
  ylab = 'Errors [a.u.]',
  plotBA = FALSE,
  topMar = 3,
  nclass = 33,
  ylim = c(-0.05,0.05),
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
png(file = paste0(figDir,'/distZscoresNew.png'),
    width= gPars$reso,height=gPars$reso)
ErrViewLib::plotDistHist(
  calc, z,
  main = '',
  plotGauss = TRUE,
  plotReg = FALSE, degree = 1,
  xlab = 'Calculated value [a.u.]',
  ylab = 'z-scores',ylim = c(-5,5),
  plotBA = FALSE,
  topMar = 3,
  nclass = 33,
  gPars  = gPars,
  scalePoints = 0.15
)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()

label = 3
png(file = paste0(figDir,'/LZV1_Example.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  calc,
  z,
  nBin = 5, #slide = TRUE,
  ylim = c(0.5,1.5),
  xlab  = 'Calculated value [a.u.]',
  label = label,
  gPars = gPars
)
dev.off()

label = 4
png(file = paste0(figDir,'/LZV2_Example.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  uE,
  z,
  logX = FALSE,
  nBin = 5, #slide = TRUE,
  ylim = c(0.5,1.5),
  xlab  = 'Prediction uncertainty [a.u.]',
  label = label,
  gPars = gPars
)
dev.off()

