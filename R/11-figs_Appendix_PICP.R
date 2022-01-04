if(!exists('doneSetup'))
  source('00-setup.R')


png(file = paste0(figDir,'/PICP_CI.png'),
    width = 2*gPars$reso, height = 2*gPars$reso)
par(
  mfrow=c(2,2),
  pty = gPars$pty,
  mar = c(3,3,1.5,1),
  mgp = gPars$mgp,
  tcl = gPars$tcl,
  lwd = gPars$lwd,
  cex = 1.25 *gPars$cex
)

Nseq = seq(20,300,by=5)
pseq = c(0.91,0.95,0.99,1.0)

methods = c("wilson", "wilsoncc",
            "agresti-coull",
            "clopper-pearson",
            "jeffreys")
nm = length(methods)

label = 0
for (kp in seq_along(pseq)) {
  nu = pseq[kp]

  for(j in seq_along(methods)) {
    meth = methods[j]

    res0 = matrix(NA,ncol=3,nrow=length(Nseq))
    for(i in seq_along(Nseq))
      res0[i,] = DescTools::BinomCI(nu*Nseq[i],Nseq[i],method = meth)

    matplot(Nseq, res0[,2:3],
            type = 'l', lty=1, col=j, lwd = 2* gPars$lwd,
            xlab = 'M', xaxs='i',
            ylab = '95% CI limits', ylim = c(0.8,1), yaxs = 'i',
            add = j>1,
            main = bquote(nu == .(nu)),
            cex.main=1
    )
  }
  grid()
  abline(h=0.95, lty=2)
  if(kp == 1)
  legend('left',
         bty ='n',
         legend = c('W','Wcc','AC','CP','J'),
         cex = 1,
         col = 1:nm,
         lty=1,
         lwd = 2*gPars$lwd)
  box()
  label = label + 1
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)

}
dev.off()

# Coverage and FP ####

nMC = 1e6
p0  = 0.95
CImeths = c('wilsoncc','clopper-pearson')
Nseq = seq(20,400,by=1)
res = matrix(NA,nrow=length(Nseq),ncol=1+length(CImeths))
for(j in seq_along(Nseq)) {
  N = Nseq[j]
  test  = rbinom(nMC,N,p0) / N
  for(k in seq_along(CImeths)) {
    ci = DescTools::BinomCI(p0*N, N, method = CImeths[k])
    res[j,k] =  sum(test >= ci[2] & test <= ci[3])/nMC
  }
}

pseq1 = c(0.85, 0.91, 0.99)
res2 = matrix(NA,ncol = length(pseq), nrow = length(Nseq))
for(i in seq_along(pseq1)) {
  res2[,i] = EnvStats::propTestPower(
    n.or.n1  = Nseq,
    p.or.p1  = pseq1[i],
    p0.or.p2 = 0.95,
    approx   = FALSE # use Clopper-Pearson "exact"
  )$power
}

pseq2 = c( 0.55, 0.6, 0.7)
res3 = matrix(NA,ncol = length(pseq), nrow = length(Nseq))
for(i in seq_along(pseq)) {
  res3[,i] = EnvStats::propTestPower(
    n.or.n1  = Nseq,
    p.or.p1  = pseq2[i],
    p0.or.p2 = 0.5,
    approx   = FALSE # use Clopper-Pearson "exact"
  )$power
}



png(file = paste0(figDir,'/BPCI_compare.png'),
    width = 3*gPars$reso, height = gPars$reso)
par(
  mfrow=c(1,3),
  pty = gPars$pty,
  mar = c(3,3,1.5,1),
  mgp = gPars$mgp,
  tcl = gPars$tcl,
  lwd = gPars$lwd,
  cex = 1.25*gPars$cex
)

label=1
matplot(Nseq, res, type = 'l',
        lty = 1, pch=0:3,
        col = gPars$cols[c(2,6)],
        lwd = 2*gPars$lwd,,
        xlab = 'M',
        ylim = c(0.85,0.975),
        ylab = 'Effective coverage of 95 % CI')
grid()
abline(h=0.95,lty=2)
legend('bottomright', bty ='n',
       legend = c('Wcc','CP'),
       lty = 1,
       lwd = 2*gPars$lwd,
       col = gPars$cols[c(2,6)],
       pch = NA)
box()
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.3)



label=2
matplot(Nseq,res3, type = 'l',
        lty = 1,
        col = gPars$cols[c(1,3,6)],
        lwd = 1.5*gPars$lwd,,
        xlab = 'M',
        ylim = c(0,1), ylab = 'Power', yaxs = 'i',
        main = 'p = 0.5')
grid()
abline(h=0.8,lty=2)
legend('bottomright', bty ='n',
       title = expression(nu),
       legend = pseq2,
       lty = 1,
       lwd = 2*gPars$lwd,
       col = gPars$cols[c(1,3,6)],
       pch = NA)
box()
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.3)

label=3
matplot(Nseq,res2, type = 'l',
        lty = 1,
        col = gPars$cols[c(1,3,6)],
        lwd = 1.5*gPars$lwd,,
        xlab = 'M',
        ylim = c(0,1), ylab = 'Power', yaxs = 'i',
        main = 'p = 0.95')
grid()
abline(h=0.8,lty=2)
legend('bottomright', bty ='n',
       title = expression(nu),
       legend = pseq1,
       lty = 1,
       lwd = 2*gPars$lwd,
       col = gPars$cols[c(1,3,6)],
       pch = NA)
box()
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.3)

dev.off()

