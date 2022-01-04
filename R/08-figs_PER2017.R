if(!exists('doneSetup'))
  source('00-setup.R')

# Data ----
D = read.csv('../data/PER2017_Data.csv', header = FALSE,
             check.names = FALSE, stringsAsFactors = FALSE)
systems = D[,1]
data = D[,2]
ref = D[,3]
err = ref-data
methods = colnames(err)

D = data
R = ref
E = err
N = length(R)

df   = as.data.frame(cbind(D,E))
up0 = sd(E)

regLct = lm(E ~ 0+D,data = df)
ELct = residuals(regLct)
upLct = sd(ELct)
z = ELct/upLct


xlab = expression(Calc.~frequency~group("[",cm^{-1},"]"))


label = 1
png(file = paste0(figDir,'/dist_Scale_PER2017.png'),
    width = gPars$reso, height = gPars$reso)
ErrViewLib::plotDistHist(
  D, ELct,
  plotGauss = TRUE,
  plotReg = TRUE, degree = 1,
  plotBA  = TRUE,
  ylim = range(ELct),
  xlab = xlab,
  ylab = expression(Error~group("[",cm^{-1},"]")),
  # colPoints = gPars$cols_tr2[cols],
  scalePoints = 0.2,
  topMar = 3,
  nclass=33,
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
png(file = paste0(figDir,'/LZV_PER2017.png'),
    width = gPars$reso, height = gPars$reso)
res = ErrViewLib::plotLZV(
  D, z,
  method = 'cho',
  xlab = xlab,
  label     = label,
  gPars     = gPars
)
dev.off()

label = 3
png(file = paste0(figDir,'/LCP_EQ_PER2017.png'),
    width = gPars$reso, height = gPars$reso)
Data = as.data.frame(cbind(D,R=ref))
res = ErrViewLib::plotPcoverage(
  Data,
  corTrend = TRUE,
  fo = E~ 0 + D,
  CImeth = 'eq',
  # title = paste0('Scale/EQ'),
  xlab = xlab,
  mycols    = c(2,5,7),
  ylim=c(0.2,1),
  label = label,
  gPars=gPars)
dev.off()

