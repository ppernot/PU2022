if(!exists('doneSetup'))
  source('00-setup.R')

gPars$cex = 5.5
gPars$lwd = 7

N = 10^4

C  = runif(N,-1,1)
z0 = rnorm(N)
E0 = rnorm(N)
few = sample.int(N,100)

for(m in c(-0.5, 0, 0.5)) {
  for(s in c(0.5, 1, 2)) {
    E  = rnorm(N) + m/s
    R  = E + C
    uE = rep(1/s,N)
    z  = E / uE
    tag = paste0('m=',m,'_s=',s)
    # gPars$mar = c(3,3,1.6,0.5)

    png(file = paste0(figDir,'/ATLAS_ecdf_',tag,'.png'),
        width = gPars$reso, height = gPars$reso)
    title = paste0('z~N(',m,',',s,')')
    ErrViewLib::plotUncEcdf(
      cbind(z,z0), xmin = -5, xmax = 5,
      xlab = 'z-score',
      absErrors = FALSE,
      show.Q95 = FALSE, show.leg = FALSE,
      col.index = c(5,2),
      gPars = gPars)
    legend('topleft', cex=0.75, bty='n',
           legend = c('Stand. Norm.','Sample'),
           lwd = 8, col = gPars$cols[c(2,5)]
    )
    # title(main = title)
    dev.off()
    png(file = paste0(figDir,'/ATLAS_pit_',tag,'.png'),
        width = gPars$reso, height = gPars$reso)
    ErrViewLib::plotPIT(
      z, gPars = gPars)
    dev.off()
    png(file = paste0(figDir,'/ATLAS_pp_',tag,'.png'),
        width = gPars$reso, height = gPars$reso)
    ErrViewLib::plotPpnorm(
      z, plotCI = FALSE, score = FALSE, gPars  = gPars)
    dev.off()
    png(file = paste0(figDir,'/ATLAS_pcov_',tag,'.png'),
        width = gPars$reso, height = gPars$reso)
    # gPars$mar = c(3,3,1.6,2)
    ErrViewLib::plotPcoverage(
      data.frame(R=R, D=C, uP=uE),
      xlab = 'Calculated value [a.u.]',
      mycols    = c(2,5,7),
      ylim      = c(0.,1.0),
      gPars     = gPars
    )
    dev.off()
  }
}
