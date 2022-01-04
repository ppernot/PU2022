if(!exists('doneSetup'))
  source('00-setup.R')

D = read.csv('../data/Data_wo1.csv',
             check.names = FALSE, stringsAsFactors = FALSE)
systems = D[,1]
data = D[,c(-1,-2)]
ref = D[,2]
err = ref-data
methods = colnames(err)


# Stats table ----
set.seed(123)
nBoot = 1000
M = nrow(err)

msetab   = umsetab   = pmsetab   = list()
I95lotab = uI95lotab = pI95lotab = list()
I95uptab = uI95uptab = pI95uptab = list()
u95tab   = uu95tab   = pu95tab   = list()
uptab    = uuptab    = puptab    = list()
ktab     = uktab     = pktab     = list()

for (meth in methods) {
  pmeth = sub('ccpvdz_','',meth)
  pmeth = sub('_wo','',pmeth)
  X = err[,meth]

  s1 = s2 = s3 = s4 = vector(length = nBoot)
  for(iB in 1:nBoot) {
    Xb = sample(X, M, replace = TRUE)
    s1[iB] = ErrViewLib::mse(Xb)
    s2[iB] = ErrViewLib::hd(Xb, q = 0.025)
    s3[iB] = ErrViewLib::hd(Xb, q = 0.975)
    s4[iB] = ErrViewLib::rmsd(Xb)
  }
  msetab[[pmeth]]  = mean(X)
  umsetab[[pmeth]] = sd(s1)
  pmsetab[[pmeth]] = ErrViewLib::prettyUnc(
    msetab[[pmeth]],umsetab[[pmeth]],numDig = 1)

  I95lotab[[pmeth]]  = ErrViewLib::hd(X, q = 0.025)
  uI95lotab[[pmeth]] = sd(s2)
  pI95lotab[[pmeth]] = ErrViewLib::prettyUnc(
    I95lotab[[pmeth]],uI95lotab[[pmeth]],numDig = 1)

  I95uptab[[pmeth]]  = ErrViewLib::hd(X, q = 0.975)
  uI95uptab[[pmeth]] = sd(s3)
  pI95uptab[[pmeth]] = ErrViewLib::prettyUnc(
    I95uptab[[pmeth]],uI95uptab[[pmeth]],numDig = 1)

  u95tab[[pmeth]] = 0.5*(I95uptab[[pmeth]]-I95lotab[[pmeth]])
  uu95tab[[pmeth]] = 0.5*sd(s3-s2)
  pu95tab[[pmeth]] = ErrViewLib::prettyUnc(
    u95tab[[pmeth]],uu95tab[[pmeth]],numDig = 1)

  uptab[[pmeth]] = ErrViewLib::rmsd(X)
  uuptab[[pmeth]] = sd(s4)
  puptab[[pmeth]] = ErrViewLib::prettyUnc(
    uptab[[pmeth]],uuptab[[pmeth]],numDig = 1)

  ktab[[pmeth]] = u95tab[[pmeth]]/uptab[[pmeth]]
  uktab[[pmeth]] = sd(0.5*(s3-s2)/s4)
  pktab[[pmeth]] = ErrViewLib::prettyUnc(
    ktab[[pmeth]],uktab[[pmeth]],numDig = 1)

}
df = data.frame(
  MSE    = unlist(pmsetab),
  I95_lo = unlist(pI95lotab),
  I95_up = unlist(pI95uptab),
  u_95   = unlist(pu95tab),
  u_p    = unlist(puptab),
  k      = unlist(pktab)
)
knitr::kable(df)

sink(paste0(tabDir,'calib_DML.tex'))
print(
  xtable::xtable(
    df,
    type = 'latex',
    caption = 'Calibration statistics for case DML',
    label = "tab:cal_DML"
  ),
  comment = FALSE,
  caption.placement ='bottom'
)
sink()

