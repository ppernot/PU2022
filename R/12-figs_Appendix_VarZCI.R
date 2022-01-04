if(!exists('doneSetup'))
  source('00-setup.R')

# Compare la couverture des intervalles à 95%
library(doParallel)
registerDoParallel(detectCores()-1)
library(foreach)

genZ = function(N, pdf=c("norm","pexp","stud","unif"), nu = 2) {
  switch (
    pdf,
    norm = rnorm(N),
    unif = (runif(N)-0.5) * 2*sqrt(3),
    stud = rt(N, df=nu) * sqrt(nu-2)/sqrt(nu),
    pexp = normalp::rnormp(N, p = nu) /
      sqrt(nu^(2/nu)*gamma(3/nu)/gamma(1/nu))
  )
}
varvar = function (Z) {
  # Use formula (1) in CHo2008
  N = NROW(Z)
  mu = moments::all.moments(Z ,order.max=4, central=TRUE)[2:5]
  Vcho = (mu[4] - (N-3)/(N-1) * mu[2]^2 ) / N
  return(Vcho)
}


# Var testing ####


fileRes = paste0(tabDir,'saveRes_BasicBCaPerc.csv')

if(!file.exists(fileRes)) {

  # Generate results ###

  nMC   = 1000
  nBoot = 1500


  pdfs = c('pexp','pexp','norm','stud','unif')[1:3]
  shap = c(     1,     4,    NA,     4,    NA)[1:3]


  Nseq = c(
    seq(50, 100, by = 10),
    seq(150, 500, by = 50),
    750, 1000)

  Nseq = c(500,750,1000)

  res = data.frame(
    pdf = NA,
    nu = NA,
    N = NA,
    V = NA,
    VV = NA,
    loBBbasic = NA,
    upBBbasic = NA,
    loBBbca = NA,
    upBBbca = NA,
    loBBperc = NA,
    upBBperc = NA
  )

  for (iN in seq_along(Nseq)) {
    N = Nseq[iN]

    for (iS in seq_along(pdfs)) {
      pdf = pdfs[iS]
      nu  = shap[iS]

      print(c(N, pdf, nu))

      resMC = foreach(
        i = 1:nMC,
        .combine = rbind,
        .packages = c('Hmisc', 'boot')
      ) %dopar% {
        Z     = genZ(N, pdf, nu)
        bst = boot::boot(Z, Hmisc::wtd.var, R= nBoot, stype = "f",
                         parallel = 'no', normwt = TRUE)
        bci = boot::boot.ci(bst, conf=0.95, type = c("bca","perc","basic") )
        data.frame(
          pdf,
          nu ,
          N,
          V = var(Z),
          VV = varvar(Z),
          loBBbasic = bci$basic[1,4],
          upBBbasic = bci$basic[1,5],
          loBBbca   = bci$bca[1,4],
          upBBbca   = bci$bca[1,5],
          loBBperc  = bci$perc[1,4],
          upBBperc  = bci$perc[1,5]
        )
      }
      res = rbind(res, resMC)
      if(iN == 1 & iS == 1)
        res = res[-1, ] # Remove dummy first line
      write.csv(res,
                file = fileRes,
                row.names = FALSE)
    }
  }
}


# Plots ####


resp = read.csv(file = fileRes)

Nseq = sort(unique(resp$N))

for(iS in seq_along(pdfs)) {
  pdf = pdfs[iS]
  nu  = shap[iS]

  Vm = VVm =
    loBBmbasic = upBBmbasic = pBBbasic =
    loBBmbca   = upBBmbca   = pBBbca   =
    loBBmperc  = upBBmperc  = pBBperc  =
    pChisq = pCho = pEmp = loEmp = upEmp =
    loChom = upChom = loChisq= upChisq = c()
  for(iN in seq_along(Nseq)) {
    N = Nseq[iN]

    # Extract data corresponding to pdf
    if(!is.na(nu))
      sel = resp$pdf == pdf & resp$nu == nu & resp$N == N
    else
      sel = resp$pdf == pdf & resp$N == N
    subres = resp[sel,]

    V    = as.numeric(subres$V)
    VV   = as.numeric(subres$VV)
    loBBbasic = as.numeric(subres$loBBbasic)
    upBBbasic = as.numeric(subres$upBBbasic)
    loBBbca   = as.numeric(subres$loBBbca)
    upBBbca   = as.numeric(subres$upBBbca)
    loBBperc  = as.numeric(subres$loBBperc)
    upBBperc  = as.numeric(subres$upBBperc)

    # Mean values and PICP
    ## Normal approx.
    Lo = qchisq(0.025, df=N-1) / (N-1)
    Up = qchisq(0.975, df=N-1) / (N-1)
    loChisq[iN] = Lo
    upChisq[iN] = Up
    pChisq[iN]  = mean(V >= Lo & V <= Up); print(c(N,pdf,nu,pChisq[iN]))

    ## Cho2008
    SD = sqrt(VV)
    limsN      = cbind(-1.96 * SD, 1.96 * SD) + V
    Vm[iN]     = mean(V)
    VVm[iN]    = mean(VV)
    loChom[iN] = mean(limsN[,1])
    upChom[iN] = mean(limsN[,2])
    pCho[iN]   = mean(limsN[,1] <= 1 & limsN[,2] >= 1)

    ## Bootstrap
    loBBmbasic[iN] = mean(loBBbasic)
    upBBmbasic[iN] = mean(upBBbasic)
    pBBbasic[iN]   = mean(loBBbasic <= 1 & upBBbasic >= 1)

    loBBmbca[iN] = mean(loBBbca)
    upBBmbca[iN] = mean(upBBbca)
    pBBbca[iN]   = mean(loBBbca <= 1 & upBBbca >= 1)

    loBBmperc[iN] = mean(loBBperc)
    upBBmperc[iN] = mean(upBBperc)
    pBBperc[iN]   = mean(loBBperc <= 1 & upBBperc >= 1)

    ## Empirical value
    qEmp = ErrViewLib::vhd(V)
    loEmp[iN] = qEmp[1]
    upEmp[iN] = qEmp[2]
    pEmp[iN]  = mean( qEmp[1] <= V & qEmp[2] >= V)

  }

  tag = paste0(
    # ciMeth,'_',
    pdf,'_',
    ifelse(is.na(nu),'',nu)
  )
  png(file = paste0(figDir,'/zVarCI_',tag,'.png'),
      width = 1.75*gPars$reso, height = gPars$reso)
  par(mfrow = c(1,2),
      mar = gPars$mar,
      mgp = gPars$mgp,
      pty = 's',
      tcl = gPars$tcl,
      cex = gPars$cex,
      lwd = 2*gPars$lwd)
  matplot(
    Nseq,
    cbind(
      # Vm,
      loChisq, upChisq,
      loChom,  upChom,
      loBBmbasic, upBBmbasic,
      loBBmbca  , upBBmbca,
      loBBmperc , upBBmperc,
      loEmp,   upEmp
    ),
    type = 'l',
    lty = 1, lwd = 2*gPars$lwd,
    col = gPars$cols[c(2,2,3,3,4,4,5,5,6,6,7,7)],
    xlab = 'N', log ='x',
    ylab = 'z-score variance CI95 limits'
  )
  grid()
  abline(h=1,col='gray50')
  box()
  label = 1
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)

  matplot(
    Nseq,
    cbind(
      pChisq,
      pCho,
      pBBbasic,
      pBBbca,
      pBBperc,
      pEmp
    ),
    type = 'l',
    lty = 1, lwd = 2*gPars$lwd,
    col = gPars$cols[2:7],
    xlab = 'N', log ='x',
    ylab = 'PICP', yaxs = 'i',
    ylim = c(0.75,1)
  )
  grid()
  abline(h=0.95,col = 'gray50')
  legend(
    'bottomright', bty = 'y', bg = 'white', box.col = 'white',
    legend = c('Chisq','Cho','Basic','BCa','Perc','Emp.'),
    lty = 1, pch = NA,
    col = gPars$cols[2:7]
  )
  box()
  label = 2
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = gPars$cex,
    line = 0.3)
  dev.off()
}

