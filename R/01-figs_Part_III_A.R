if(!exists('doneSetup'))
  source('00-setup.R')


set.seed(1)

N = 1000
ref  = rnorm(N)
err0 = rnorm(N,0,0.1)

# Biased data
calc1 = 1.25 * ref + err0
err1  = ref - calc1
uErr1 = sd(err1)
uC1 = rep(uErr1, length(calc1))
z1 = err1 / uErr1

# Unbiased data
reg = lm(ref ~ calc1)
calc2 = predict(reg)
err2 =  residuals(reg)
uErr2 = sd(err2)
uC2 = rep(uErr2, length(calc2))
z2 = err2 / uErr2


# Fig example ----

label = 1
png(file = paste0(figDir,'/biasedErrors.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotDistHist(
  calc1, err1,
  main = '',
  plotGauss = TRUE,
  plotReg = FALSE, degree = 1,
  xlab = 'Calculated value',
  ylab = 'Error',
  plotBA = TRUE,
  topMar = 3,
  gPars = gPars,
  scalePoints = 0.1
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
png(file = paste0(figDir,'/PIT_biasedErrors.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotPIT(
  z1,
  title = '',
  label = label,
  gPars = gPars
)
dev.off()

label = 3
png(file = paste0(figDir,'/ppNorm_biasedErrors.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotPpnorm(
  z1,
  title = '',
  plotCI = TRUE,
  label = label,
  gPars = gPars
)
dev.off()

selPos = calc1>=0

label = 4
png(file = paste0(figDir,'/biasedErrorsPartial.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotDistHist(
  calc1, err1,
  main = '',
  plotGauss = TRUE,
  plotReg = FALSE, degree = 1,
  xlab = 'Calculated value',
  ylab = 'Error',
  plotBA = TRUE,
  topMar = 3,
  select = selPos,
  gPars  = gPars,
  scalePoints = 0.1
)
if(label > 0)
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
dev.off()

label = 5
png(file = paste0(figDir,'/PIT_biasedErrorsPartial.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotPIT(
  z1[selPos],
  title = '',
  label = label,
  gPars = gPars
)
dev.off()

label = 6
png(file = paste0(figDir,'/ppNorm_biasedErrorsPartial.png'),
    width=gPars$reso,height=gPars$reso)
ErrViewLib::plotPpnorm(
  z1[selPos],
  title = '',
  plotCI = TRUE,
  label = label,
  gPars = gPars
)
dev.off()

label = 1
png(file = paste0(figDir,'/pCoverage_biasedErrors.png'),
    width = gPars$reso, height = gPars$reso)
Data = as.data.frame(cbind(D=calc1, R=ref, uP=uC1))
ErrViewLib::plotPcoverage(
  Data,
  mycols    = c(7, 5, 2),
  ylim      = c(0.,1.0),
  title     = '',
  label     = label,
  gPars     = gPars
)
dev.off()

label = label +1
png(
  file = paste0(figDir, '/pCoverage_unbiasedErrors.png'),
  width = gPars$reso,
  height = gPars$reso
)
Data = as.data.frame(cbind(D = calc2, R = ref, uP = uC2))
ErrViewLib::plotPcoverage(
  Data,
  mycols    = c(7, 5, 2),
  ylim      = c(0., 1.0),
  label     = label,
  gPars     = gPars
)
dev.off()
