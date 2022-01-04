if(!exists('doneSetup'))
  source('00-setup.R')

# BEEF Calibration ----
fileName = '../data/PAN2015_Data.csv'
table = read.table(
  fileName,
  header = FALSE,
  stringsAsFactors = FALSE,
  encoding = 'utf8',
  colClasses = c(
    'character',
    'numeric',
    'numeric',
    'numeric',
    'numeric',
    'numeric',
    'character',
    'numeric',
    'numeric',
    'numeric',
    'numeric',
    'numeric'
  )
)
# reshape table 12 columns -> 6 columns
tab = table[, 1:6]
colnames(tab) = paste0('V', 7:12)
tab = rbind(tab, table[, 7:12])
# Get rid of last line (duplicated to fill original table)
tab = tab[-nrow(tab), 1:4]
colnames(tab) = c('Species', 'Expt', 'BEEF', 'uBEEF')
ord = order(tab[, 3])
Calc  = tab[ord, 3]
uCalc = tab[ord, 4]
Ref   = tab[ord, 2]
err   = Ref - Calc
z     = err / uCalc
N     = length(z)

reg = lm(Ref~Calc)
lc_err = residuals(reg)
uErr = sd(lc_err)
lc_calc = predict(reg)
lc_z  = lc_err / uCalc
lc_z0 = lc_err / uErr

# Check consistency of z-scores
print(ErrViewLib::varZCI(z))


label = 1
png(file = paste0(figDir,'/dist_PAN2015.png'),
    width =gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  Calc,err,uCalc,
  plotGauss = TRUE,
  plotReg = TRUE, degree = 1,
  ylim = range(c(err-2*uCalc,err+2*uCalc)),
  xlab = 'Calculated heat [kcal/mol]',
  ylab = 'Error [kcal/mol]',
  nclass = 17, plotBA = TRUE,
  topMar = 3,
  scalePoints = 0.25,
  gPars=gPars)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()

label = 2
png(file = paste0(figDir,'/dist_zscores_PAN2015.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  Calc,z,
  plotGauss = TRUE,
  plotReg = TRUE, degree = 1,
  xlab = 'Calculated heat [kcal/mol]',
  ylab = 'z-score',
  nclass= 17,plotBA = TRUE,
  topMar = 3,
  scalePoints = 0.25,
  gPars=gPars)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()

label = 3
png(file = paste0(figDir,'/LZV1_PAN2015.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  Calc,
  z,
  logX  = FALSE,
  nBin  = 5,
  slide = TRUE,
  xlab  = 'Calculated heat [kcal/mol]',
  label = label,
  gPars = gPars
  )
dev.off()

label = 4
png(file = paste0(figDir,'/LZV2_PAN2015.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotLZV(
  uCalc,
  z,
  # logX = TRUE,
  nBin = 5,
  slide = TRUE,
  xlab  = 'Prediction uncertainty [kcal/mol]',
  label = label,
  gPars = gPars
)
dev.off()


# A-posteriori analysis ####

label = 1
png(file = paste0(figDir,'/dist_PCT1_PAN2015.png'),
    width =gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  Calc,lc_err,uErr,
  plotGauss = TRUE,
  # plotReg = TRUE, degree = 1,
  ylim = range(c(lc_err-2*uErr,lc_err+2*uErr)),
  xlab = 'Calculated heat [kcal/mol]',
  ylab = 'Error [kcal/mol]',
  nclass = 17, plotBA = TRUE,
  topMar = 3,
  scalePoints = 0.25,
  gPars=gPars)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()


label = label + 1
png(file = paste0(figDir,'/pCoverage_LCT1_PRED_PAN2015.png'),
    width = gPars$reso, height = gPars$reso)
res = ErrViewLib::plotPcoverage(
  as.data.frame(cbind(D=Calc, R=Ref)),
  mycols    = c(2,5,7),
  corTrend  = TRUE,
  fo        = E ~ D,
  CImeth    = 'pred',
  ylim      = c(0.4,1.0),
  # title     = 'PTC1/PRED calibration',
  xlab      = 'Calculated heat [kcal/mol]',
  label     = label,
  gPars     = gPars
)
dev.off()

label = label + 1
png(file = paste0(figDir,'/pCoverage_LCT1_EQ_PAN2015.png'),
    width = gPars$reso, height = gPars$reso)
res = ErrViewLib::plotPcoverage(
  as.data.frame(cbind(D=Calc, R=Ref)),
  mycols    = c(2,5,7),
  corTrend  = TRUE,
  fo        = E ~ D,
  CImeth    = 'eq',
  ylim      = c(0.4,1.0),
  # title     = 'PTC1/PRED calibration',
  xlab      = 'Calculated heat [kcal/mol]',
  label     = label,
  gPars     = gPars
)
dev.off()






