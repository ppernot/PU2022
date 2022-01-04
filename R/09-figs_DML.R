if(!exists('doneSetup'))
  source('00-setup.R')

# Delta-ML ----

D = read.csv('../data/Data_wo1.csv',
             check.names = FALSE, stringsAsFactors = FALSE)
systems = D[,1]
data = D[,c(-1,-2)]
ref = D[,2]
err = ref-data
methods = colnames(err)
R = ref
N = length(ref)

niceMethNames = c(
  'HF','MP2','SLATM-L2',
  expression(Delta-ML/SLATM-L2),
  expression(Delta-ML/CM-L1)
)

# PTC1/PRED
label = 0
for (iMeth in c(2,3,5)) {
  label = label + 1

  set.seed(123)
  meth = methods[iMeth]
  pmeth = sub('_wo','',sub('ccpvdz_','',meth))
  methScore = 'PRED'

  z = zEst(
    Data = as.data.frame(cbind(D=data[,meth],R=R)),
    corTrend = TRUE,
    fo = E~D,
    UQmeth = 'pred'
  )$z

  png(file = paste0(figDir,'/CalCurve_',pmeth,'_',methScore,'.png'),
      width = gPars$reso, height = gPars$reso)
  ErrViewLib::plotPpnorm(
    z,
    title = niceMethNames[iMeth],
    gPars = gPars,
    label = label)
  dev.off()

  C = data[,meth]
  png(file = paste0(figDir,'/LCP_',pmeth,'_',methScore,'.png'),
      width = gPars$reso, height = gPars$reso)
  ErrViewLib::plotPcoverage(
    as.data.frame(cbind(R=R, D=C)),
    CImeth    = 'pred',
    corTrend  = TRUE,
    fo        = E~D,
    mycols    = c(2,5,7),
    xlim      = range(C),
    ylim      = c(0.3,1.0),
    title     = niceMethNames[iMeth],
    xlab      = 'Calculated EAE [kcal/mol]',
    label     = label,
    gPars     = gPars
  )
  dev.off()

  png(file = paste0(figDir,'/LZV_',pmeth,'_',methScore,'.png'),
      width = gPars$reso, height = gPars$reso)
  ErrViewLib::plotLZV(
    data[,meth],
    z,
    method = 'cho',
    slide = TRUE,
    xlab  = 'Calculated EAE [kcal/mol]',
    title = niceMethNames[iMeth],
    label = label,
    gPars = gPars
  )
  dev.off()
}
