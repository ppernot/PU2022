# lib = "ErrViewLib"
# if(!require(lib,character.only = TRUE))
# devtools::install_github(paste0("ppernot/",lib))
# library(lib,character.only = TRUE)
#
library(CHNOSZ)

figDir = '../results/figs/'
tabDir = '../results/tables/'

gPars = ErrViewLib::setgPars(type='publish')
gParsScreen = ErrViewLib::setgPars()

BPCImethod = 'wilsoncc'

# Functions ----

# frms = function(x,indx=1:length(x),...) {
#   sqrt(mean(x[indx]^2))
# }
# fp95 = function(x, indx=1:length(x),eq95,...) {
#   mean( x[indx] >= eq95[1] & x[indx] <= eq95[2])
# }
# fqrat = function(x,indx=1:length(x),...) {
#   diff(vhd(x[indx])) / (2*1.96*sd(x[indx]))
# }
# trapz = function (
#   x,
#   y
# ) {
#   if(length(x) == 0)
#     return(0)
#   idx = 2:length(x)
#   return(
#     as.double((x[idx] - x[idx - 1]) %*%
#                 (y[idx] + y[idx - 1]))/2
#   )
# }
#
# plotReliability = function(
#   E, uC,
#   nBoot = 1000,
#   nBin = 10,
#   title = '',
#   gPars
# ) {
#
#   for (n in names(gPars))
#     assign(n, rlist::list.extract(gPars, n))
#
#   par(
#     mfrow = c(1, 1),
#     mar = mar,
#     mgp = mgp,
#     pty = 's',
#     tcl = tcl,
#     cex = cex,
#     lwd = lwd
#   )
#
#   ord = order(uC)
#   uOrd = uC[ord]
#   eOrd = E[ord]
#
#   p = seq(0,1,length.out=nBin+1)[-1]
#   br = quantile(uOrd,p=p)
#   cl = c()
#   for(i in seq_along(uOrd))
#     cl[i] = which(br >= uOrd[i])[1]
#
#   rmv = rmse = maz = urmv = urmse = c()
#   for(i in seq_along(br)){
#     sel = cl==i
#     rmv[i]  = frms(uOrd[sel])
#     urmv[i] = sd(boot::boot(uOrd[sel], frms, R=nBoot)$t)
#     rmse[i] = frms(eOrd[sel])
#     urmse[i]= sd(boot::boot(eOrd[sel], frms, R=nBoot)$t)
#   }
#   rmv0  = sqrt(mean(uOrd^2))
#   rmse0 = sqrt(mean(eOrd^2))
#
#   xlim = ylim = range(c(rmv, rmse))
#   plot( rmse, rmv,
#         xlab = 'RMSE', xlim = xlim,
#         ylab = 'RMV',  ylim = ylim,
#         pch = 19, col = cols[5],
#         main = title)
#   grid()
#   segments(rmse, rmv - 2 * urmv, rmse, rmv + 2 * urmv,   col = cols[5])
#   segments(rmse - 2 * urmse, rmv, rmse + 2 * urmse, rmv, col = cols[5])
#   abline(a = 0,b = 1, lty = 2, col = cols[1])
#   abline(h = rmv0, lty = 3, col = cols[3])
#   abline(v = rmse0, lty = 3, col = cols[3])
#   box()
#
# }
#
# zLinReg = function(
#   Ref, Calc,
#   deg = 0
# ){
#
#   # Predictor
#   x = Calc
#   # Errors
#   y = Ref - Calc
#
#   # Linear Trend correction
#   fo = y ~ 1
#   if (deg > 0)
#     fo = as.formula(
#       paste0('y ~ 1 +',
#              paste0(
#                'I(x^', 1:deg, ')',
#                collapse = '+'
#              )))
#   reg = lm(fo)
#   E   = residuals(reg)
#
#   # Prediction uncertainty
#   pE  = predict(reg, se.fit = TRUE)
#   # Ensure that the prediction variance is
#   # equal to the residuals variance
#   uE  = sqrt( pE$se.fit^2 + var(E) - mean(pE$se.fit^2))
#
#   # z-scores
#   z   = E / uE
#
#   return(
#     list(
#       deg = deg,
#       E   = E,
#       uE  = uE,
#       z   = z
#     )
#   )
# }
#
zEst = function(
  Data,
  corTrend = FALSE,
  fo       = NA,
  UQmeth   = c('seq','pred','dist')
  # meth = c('C0','LCT0','LCT1','LCT2','SEQ')
){

  meth = match.arg(UQmeth)

  E = Data$R - Data$D

  if(corTrend) {
    environment(fo) <- environment()
    reg = lm(fo, data = Data)
    E  = residuals(reg)
  }

  if(meth == 'dist') {
    E = E - mean(E)
    uE = sd(E)

  } else if (meth == 'pred') {
    # Prediction uncertainty
    pE  = predict(reg, se.fit = TRUE)$se.fit
    # Ensure that the prediction variance is
    # equal to the residuals variance
    uE  = sqrt( pE^2 + var(E) - mean(pE^2) )

  } else {
    # SEQ
    E = E - median(E)
    qE = ErrViewLib::vhd(E) # Limits of the 95% CI
    uE = diff(qE) /3.92

  }

  # z-scores
  z = E / uE

  return(
    list(
      E   = E,
      uE  = uE,
      z   = z
    )
  )
}
bootPICP = function(
  df,
  nBoot = 1000,
  numDig = 2
) {
  N = length(df$E)
  tab = rep(0,nBoot)
  for(iB in seq_along(tab)) {
    indx = sample.int(N,N,TRUE)
    tab[iB] = mean( abs(df$E[indx]) <= df$Up[indx] )
  }
  ErrViewLib::prettyUnc(mean(tab),sd(tab),numDig)
}
#
#
# calStats = function(
#   X,
#   nBoot = 0,
#   nMC   = 10^4
# ) {
#   # Calibration statistics : MisCal and CalErr
#   # with bootstrap uncertainty
#
#   M = length(X)
#
#   # Theoretical probabilities
#   q = seq(-6, 6, length.out = 100)
#   pt = c(0,pnorm(q),1)
#
#   # Monte Carlo estimation of 95% CI for normal
#   misCal = vector(length = nMC)
#   for (i in 1:nMC){
#     pe  = ecdf(rnorm(M))(q)
#     misCal[i] = trapz(pt,abs(c(0,pe,1)-pt))
#   }
#   misCalUp = ErrViewLib::hd(misCal, q = 0.975)
#
#   pe = c(0,ecdf(X)(q),1)
#   misCal = trapz(pt,abs(pe-pt))
#   calErr = sqrt(sum((pt-pe)^2))
#
#   if(nBoot > 0) {
#     s1 = s2 = vector(length = nBoot)
#     for(iB in 1:nBoot) {
#       Xb = sample(X, M, replace = TRUE)
#       peb = c(0,ecdf(Xb)(q),1)
#       s1[iB] = trapz(pt,abs(peb-pt))
#       s2[iB] = sqrt(sum((pt-peb)^2))
#     }
#     uMisCal = sd(s1)
#     uCalErr = sd(s2)
#   }
#
#   return(
#     list(
#       misCal   = misCal,
#       uMisCal  = ifelse(nBoot>0,uMisCal,NA),
#       misCalUp = misCalUp,
#       calErr   = calErr,
#       uCalErr  = ifelse(nBoot>0,uCalErr,NA)
#     )
#   )
# }
#
# picpTable = function(
#   z,
#   dist      = c('norm','t'),
#   shape     = 2,
#   prob      = c(0.50,0.75,0.95),
#   BPCImethod= 'wilsoncc'
# ) {
#
#   dist    = match.arg(dist)
#
#   # GLobal PICP
#   S = picp = upicp = lwr = upr = test = testCP = power =
#     vector(length = length(prob))
#   N = length(z)
#   for(i in seq_along(prob)) {
#
#     p = prob[i]
#     plow  = (1 - p) / 2
#     psup  = (1 + p) / 2
#     qlow  = switch(
#       dist,
#       norm = normalp::qnormp(plow, p = shape) /
#         sqrt(shape^(2/shape)*gamma(3/shape)/gamma(1/shape)),
#       t    = sqrt((shape-2)/shape) * qt(plow, df = shape)
#     )
#     qsup  = switch(
#       dist,
#       norm = normalp::qnormp(psup, p = shape) /
#         sqrt(shape^(2/shape)*gamma(3/shape)/gamma(1/shape)),
#       t    = sqrt((shape-2)/shape) * qt(psup, df = shape)
#     )
#     S[i] = sum( z >= qlow & z <= qsup )
#
#     picp[i] = S[i]/N
#     upicp[i] = sqrt(picp[i]*(1-picp[i])/N)
#
#     ci = DescTools::BinomCI(S[i], N, method = BPCImethod)
#     lwr[i] = ci[1,2]
#     upr[i] = ci[1,3]
#
#     test[i] = p >= lwr[i] & p <= upr[i]
#
#     testCP[i] = binom.test(S[i],N,p)$p.value > 0.05
#
#     power[i]  =   EnvStats::propTestPower(
#       n.or.n1  = N,
#       p.or.p1  = picp[i],
#       p0.or.p2 = p,
#       approx   = FALSE # use Clopper-Pearson "exact"
#     )$power
#   }
#
#   data.frame(
#     p = prob,
#     S = S,
#     N = N,
#     PICP = signif(picp,2),
#     uPICP = signif(upicp,2),
#     I_95_low = signif(lwr,2),
#     I_95_upr = signif(upr,2),
#     Test = test,
#     TestCP = testCP,
#     rejPower = signif(power,2)
#   )
# }

doneSetup = TRUE
