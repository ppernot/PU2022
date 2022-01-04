if(!exists('doneSetup'))
  source('00-setup.R')

# BEEF-vib ----
file='../data/PAR2019_Data.csv'
dat = read.csv(file)

method='mu'
Calc  = dat[[method]]
Ref   = dat[['expt.']]
uCalc = dat[['sigma']]
err   = Ref - Calc
z     = err / uCalc
zs = ErrViewLib::varZCI(z)
print(zs)

err_lc = err - mean(err)
z_lc   = err_lc / uCalc
uPH    = sd(err_lc)
z_PH   = err_lc / uPH

zs = ErrViewLib::varZCI(z_PH)
print(zs)

label = 1
png(file = paste0(figDir,'/dist_PAR2019.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  Calc,err,uCalc,
  plotGauss = TRUE,
  ylim = range(c(err-2*uCalc,err+2*uCalc)),
  xlab = 'Calculated frequency [meV]',
  ylab = 'Error [meV]',
  # main = 'PAR2019',
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
png(file = paste0(figDir,'/dist_zscores_PAR2019.png'),
    width = 1.1*gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  Calc,z,
  plotGauss = TRUE,
  plotReg = TRUE, degree = 1,
  xlab = 'Calculated frequency [meV]',
  ylab = 'z-score',
  plotBA = TRUE,
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
