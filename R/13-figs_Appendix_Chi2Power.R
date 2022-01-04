if(!exists('doneSetup'))
  source('00-setup.R')

library(EnvStats)

fileRes = file=paste0(tabDir,'Chi2Power.Rda')

if(!file.exists(fileRes)) {
  # Generate results ####

  power.chi2.test = function(
    N,
    mu1 = 0,
    var1 = 1,
    var0 = 1,
    reps = 25000,
    alpha = 0.05
  ) {
    sigma1 = sqrt(var1)
    pvalues = numeric(reps)
    for (i in 1:reps) {
      x1 = rnorm(N, mu1, sigma1)
      pvalues[i] = EnvStats::varTest(x1,
                                     alternative = 'two.sided',
                                     sigma.squared = var0)$p.value
    }
    mean(pvalues < alpha)
  }

  varSeq = c(seq(0.5,0.9,by=0.1),seq(1.1,1.5,by=0.1))
  nSeq = c(seq(10,90,by=10),seq(100,900,by=100),seq(1000,10000,by = 1000))

  pwr = matrix(NA,nrow=length(nSeq),ncol=length(varSeq))
  for(i in seq_along(nSeq)) {
    N = nSeq[i]
    for(j in seq_along(varSeq))
      pwr[i,j] = power.chi2.test(N, var1 = varSeq[j] )
  }

  save(varSeq,nSeq,pwr,file=fileRes)

}

# Generate fig ####

load(fileRes)

png(file = paste0(figDir,'Chi2Power.png'),
    width = gPars$reso, height = gPars$reso)
par(
  mfrow = c(1,1),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  lwd = gPars$lwd
)
lty = c(rep(1,5),rep(2,5))
col = gPars$cols[c(1:4,6,6,rev(1:4))]
matplot(
  nSeq, pwr,
  type='l',
  lty = lty, lwd = 1.75*gPars$lwd,
  col = col,
  xaxs = 'i', xlim = c(10,10^4), log = 'x', xlab = 'Sample size',
  yaxs = 'i', ylim = c(0,1), ylab = 'Power of test'
)
grid(equilogs = FALSE)
abline(h = 0.8, lty = 2, col = 'gray40', lwd = 1.5*gPars$lwd)
legend(
  'bottomright', bty = 'y',
  ncol = 1, cex = 0.75, lwd = 1.75*gPars$lwd,
  legend = paste0('Var = ', signif(varSeq,2)),
  pch= NA, lty = lty, col=col,
  bg = 'white', box.col = 'white'
)
box()
dev.off()

