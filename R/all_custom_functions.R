
# Date: 26-02-2021

#**************************#
#   all custom functions   #
#**************************#



#' These are all my custom R functions!


#*********************************************************************************
#   EMPTY CONSOLE   ####
#*********************************************************************************
lllllllll <- function(){cat(rep("\n", 50))}


#*********************************************************************************
#   SET UP GRAPHIC DEVICE   ####
#*********************************************************************************
parSave <- function(){
  parT <- par(no.readonly=TRUE)
  parT$pin <- NULL  # Do not want to save this, can only lead to problems when trying to reset to that (perhaps some other vars also)
  return(parT)
}

parReset <- function(parOld, split=FALSE){
  par(parOld)
  if(split){par(mfrow=c(2,2))}
}

parMarModify <- function(f=0.3, reduce=FALSE){
  stopifnot(0<f)
  if(reduce){
    stopifnot(f<=1)
    parM <- par()$mar
    M.d <- parM*f
    par(mar=parM-M.d)
  }else{
    parM <- par()$mar
    M.d <- parM*f
    par(mar=parM+M.d)
  }
}

setUpGraph <- function(row=2, col=2){
  oldPar <- parSave()  # In case I want to automatically save my par setup
  dev.new()
  par(mfrow=c(row,col))
  invisible(oldPar)
}


#*********************************************************************************
#   ASSIGN ARGUMENTS OF FUNCTION TO WORK INSIDE FUNCTION   ####
#*********************************************************************************
asArguments <- function(argString){
  t <- strsplit(argString, split = ', ')[[1]]   # Remove commas
  t <- t[grepl(pattern = '=', t)]   # Remove elements without "="
  eval(parse(text = t), parent.frame())   # This combination makes it possible to
                                          # perform these calls provided as strings.
}


#*********************************************************************************
#   NICE UNIVARIATE PLOT   ####
#*********************************************************************************
niceUnivPlot <- function(numVar, catVar=NULL, violin=TRUE, showMean=TRUE,
                         bw='nrd0', jitFactor=0.2, add.ylim=0, ylim.cust=NULL,
                         densScl=0.5, main=NULL, sigTest=FALSE, multCmp=FALSE){
  ### Check some requirements:
  stopifnot(is.numeric(numVar))
  stopifnot(is.null(ylim.cust) | length(ylim.cust)==2)   # Check valid entries for ylim.cust
  ### Set ylim:
  if(is.null(ylim.cust)){
    ylms <- c(min(numVar, na.rm = TRUE), max(numVar, na.rm = TRUE))
  }else{   # Custom ylim
    ylms <- ylim.cust
  }
  ylms <- ylms + c(-add.ylim, add.ylim)   # Add the add.ylim
  ### Generate catvar:
  catVar.nm <- deparse(substitute(catVar))   # Get name of catVar
  if(is.null(catVar)){
    catVar <- factor(rep(1,length(numVar)))   # 1-level factor
  }else{
    catVar <- factor(catVar)   # Dont use as.factor, otherwise non-present levels are not removed.
  }
  #**********************
  ### Run group comparisons:
  if(nlevels(catVar) > 1 & sigTest){   # Only run if catVar has entries of multiple levels
    ### List group comparisons in table:
    grInd <- data.frame(t(utils::combn(1:nlevels(catVar),2)))
    ### Apply wilcoxon tests to compare corresponding groups (collect pvalues):
    wilxP <- apply(grInd, 1, function(x){ wilcox.test(x = numVar[as.numeric(catVar)==x[1]],
                                                      y = numVar[as.numeric(catVar)==x[2]])$p.value })
    ### Combine to table (add bonferonni corrected pvals and y values in plot):
    grCmp.0 <- cbind(grInd, wilxP, 'wilxP_bonf'=wilxP*nrow(grInd))
    ### Add stars:
    sigStr <- c('', '*', '**', '***')
    sigThr <- c(1.1, 0.05, 0.01, 0.001)   # 1.1 in case pvalue is rounded to 1 (difference in next line has to be negative at least once)
    if(multCmp){
      ### Filter out only significant tests:
      grCmp <- grCmp.0[grCmp.0$wilxP_bonf < 0.05,]
      ### Have to stop in case there are no significant differences:
      if(nrow(grCmp)==0){stop('There are no significant group differences. Set sigTest to FALSE in order to draw the plot.')}
      ### Get stars according to bonferroni corrected pvalues:
      grCmp$strs <- sapply(grCmp[,'wilxP_bonf'], function(x){ sigStr[max(which((x - sigThr) < 0))] })
    }else{
      ### Filter out only significant tests:
      grCmp <- grCmp.0[grCmp.0$wilxP < 0.05,]
      ### Have to stop in case there are no significant differences:
      if(nrow(grCmp)==0){stop('There are no significant group differences. Set sigTest to FALSE in order to draw the plot.')}
      ### Get stars according to uncorrected pvalues:
      grCmp$strs <- sapply(grCmp[,'wilxP'], function(x){ sigStr[max(which((x - sigThr) < 0))] })
    }
    ### Add y-values for plot (should not be too far apart):
    grCmp$yVal <- max(ylms) + (abs(min(ylms)-max(ylms))*0.05 * 1:nrow(grCmp))
    ### Adapt ylms:
    ylms <- c(ylms[1], max(grCmp$yVal))
  }
  #**********************
  ### Get title name:
  main.nm <- deparse(substitute(main))
  ### Start plot:
  plot(x = jitter(as.numeric(catVar), factor = jitFactor),
       y = numVar,
       xlim = c(0,nlevels(catVar)+1),
       ylim = ylms,
       xaxt = 'n',
       col = catVar,
       ylab = deparse(substitute(numVar)),
       xlab = ifelse(catVar.nm=='NULL', '', catVar.nm),
       main = ifelse(main.nm=='NULL', paste0(deparse(substitute(numVar)), ' Plot'), main))
  ### Add legend:
  if(nlevels(catVar) > 1){   # Only needed with multiple levels
    legend('bottomright', legend = levels(catVar),
           pch=1,
           col=1:nlevels(catVar))
  }
  ### Add the violin lines:
  if(violin){
    L <- list()
    for(i in 1:nlevels(catVar)){
      d <- density(numVar[as.numeric(catVar)==i], na.rm = TRUE, bw = bw)
      L[[i]] <- data.frame(xd=d$x, yd=d$y)
    }
    names(L) <- levels(catVar)
    ### We have to scale the densities, need the maximum value for that:
    maxD <-  max(sapply(L, function(x){x$yd}))
    cexD <- densScl/maxD
    ### Now plot the densities:
    for(i in 1:nlevels(catVar)){
      lines(L[[i]]$yd*cexD + i, L[[i]]$xd, col=i)
      lines((-L[[i]]$yd)*cexD + i, L[[i]]$xd, col=i)
    }
  }
  ### Add the mean-value lines:
  if(showMean){
    for(i in 1:nlevels(catVar)){
      mVal <-  mean(numVar[as.numeric(catVar)==i], na.rm = TRUE)
      segments(x0 = i-0.3, y0 = mVal, x1 = i+0.3, y1 = mVal, col = i, lwd = 2)
    }
  }
  ### Draw lines of group comparisons:
  if(nlevels(catVar) > 1 & sigTest){   # Only run if comparisons were performed
    ### Draw lines:
    for(i in 1:nrow(grCmp)){
      xx <- as.numeric(grCmp[i,1:2])
      yy <- rep(grCmp[i,'yVal'], 2)
      lines(xx, yy)
      ### Find a good ticklength:
      tickL <- (abs(min(ylms)-max(ylms))*0.01)
      lines(c(xx[1], xx[1]), c(yy[1], yy[1] - tickL))
      lines(c(xx[2], xx[2]), c(yy[2], yy[2] - tickL))
      ### Add star:
      text(x = mean(xx), y = yy[1]+tickL, labels = grCmp[i,'strs'])
    }
  }
}


#*********************************************************************************
#   NICE PAIRS PLOT   ####
#*********************************************************************************
nicePairsPlot <- function(x, breaks='Sturges', density=FALSE, jitter=FALSE, jitFactor=1, loess=FALSE, swtchPan=FALSE){
  ### Panel diagonal:
  panel.diag <- function(x){
    usr <- par('usr')
    on.exit(par(usr))
    par(usr = c(usr[1], usr[2], 0, 1.5))
    h <- hist(x, breaks = breaks, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = 'cyan')
    ### Add density:
    if(density){
      tryd <- try(d <- density(x, na.rm = TRUE), silent = TRUE)
      if (class(tryd) != 'try-error') {
        d$y <- d$y/max(d$y)
        lines(d, col='cyan')
      }
    }
  }
  ### Panel top:
  panel.top <- function(x, y){
    if (jitter) {
      x <- jitter(x, factor = jitFactor)
      y <- jitter(y, factor = jitFactor)
    }
    points(x, y)
    ### Add Loess fit:
    if(loess){
      tryd <- try(lml <- suppressWarnings(loess(y ~ x, degree = 1, family = "symmetric")), silent=TRUE)
      if(class(tryd)!='try-error'){
        tempx <- data.frame(x = seq(min(x, na.rm = TRUE),
                                    max(x, na.rm = TRUE), length.out = 50))
        pred <- predict(lml, newdata = tempx)
        lines(x=tempx$x, y=pred, col='red', lty=2)
      }
    }
  }
  ### Panel bottom:
  panel.bottom <- function(x, y){
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    cort <- cor.test(x,y)
    p <- cort$p.value
    star <- symnum(p, corr = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", ""), legend = FALSE)
    txt <- paste0(round(cort$estimate, digits=2), star)
    cex <- abs(cort$estimate) * 3.5   # Set the size of numbers
    minCex <- 1   # Set the minimum size
    if(cex < minCex) {cex <- minCex}
    ### Write correlation coefficient and add stars
    text(0.5, 0.5, txt, cex = cex)
  }
  ### Plot the pairs-plot:
  if(swtchPan){
    ### Lower and upper pannel switched:
    pairs(x, diag.panel = panel.diag, lower.panel = panel.top, upper.panel = panel.bottom)
  }else{
    pairs(x, diag.panel = panel.diag, lower.panel = panel.bottom, upper.panel = panel.top)
  }
}


#*********************************************************************************
#   NICE 3D PLOT   ####
#*********************************************************************************
nice3DPlot <- function(X = NULL, whatToPlot = c('P','D','PD'), plotFit = c('no','lin','int','int2','int3'), catVar = factor(1), covMat = NULL, means = NULL, nSim = 500, pointCol = 1, pointSize = 5, pointTrans = 0.8, Dtransp_fac = 0.06, colres = 50, gridRes = 30, h = 0.3, DpointSize = 30, axesNames = NULL, axesTicks = FALSE, gridCol = 'grey', axesLeng = NULL, zoom = 1, add = FALSE){
  #*********************************************************************************
  #   TEST CONDITIONS   ####
  #*********************************************************************************
  stopifnot(whatToPlot[1] %in% c('P','D', 'PD'))
  stopifnot(plotFit[1] %in% c('no','lin','int','int2','int3'))
  stopifnot(is.factor(catVar))
  if(!all(unlist(lapply(as.data.frame(X[,1:3]), is.numeric)))){
    stop("The first three variables of X (used for plotting) have to be numeric")
  }


  #*********************************************************************************
  #   EMPTY PLOT IF NOTHING PROVIDED   ####
  #*********************************************************************************
  if(is.null(X) & is.null(means) & is.null(covMat)){
    if(is.null(axesNames)) {axesNames <- c('x','y','z')}
    if(is.null(axesLeng)) {axesLeng <- rep(1,3)}
    .coorSystem3D(axesNames, gridCol, axesLeng, axesTicks, zoom)
    warning("No data provided, so only empty coordinate system created.")
    return()
  }


  #*********************************************************************************
  #   PREPARE DATA   ####
  #*********************************************************************************
  ### If X provided:
  if (!is.null(X)){
    ### Select first three columns:
    if(ncol(X)>3){warning('More than three variables in the data. Will select the first three columns for plotting.')}
    dat <- as.data.frame(X[,1:3])
    ### Combine with other variables
    dat <- cbind(dat, "catV"=catVar, "pointC"=pointCol)
    ### Adapt pointCol if catVar supplied:
    if(nlevels(dat$catV) > 1){dat$pointC <- as.numeric(dat$catV)}
    ### Remove rows with missing values:
    noNaSel <- apply(dat, 1, function(x){all(!is.na(x))})
    dat <- dat[noNaSel,]
    if(sum(!noNaSel) > 0){warning("There are missing values in the data. Will remove the rows containing missing values.")}
    ### Generate x, y, z:
    colnames(dat)[1:3] <- c('x','y','z')
    ### Generate axesnames if not provided:
    if(is.null(axesNames)) {axesNames <- colnames(X)[1:3]}

    ### Simulate data:
  }else if(!is.null(covMat) & !is.null(means)){
    ### If no X provided:
    dat <- as.data.frame(mvtnorm::rmvnorm(n = nSim, mean = means, sigma = covMat))   # Alternative: method = 'svd'
    ### Combine with other variables
    dat <- cbind(dat, "catV"=catVar, "pointC"=pointCol)
    ### Generate x, y, z:
    colnames(dat)[1:3] <- c('x','y','z')
    ### Generate axesnames if not provided:
    if(is.null(axesNames)) {axesNames <- c('x','y','z')}
  } else {
    stop('No data (X) is provided. For a simulation either the vector of means or the covariance matrix is missing. To run a Simulation both are needed.')
  }


  #*********************************************************************************
  #   DENSITY CALCULATION   ####
  #*********************************************************************************
  ### Apply the kernel density estimator:
  dens <-  misc3d::kde3d(dat$x,dat$y,dat$z, h = h, n = gridRes)
  ### Create the matrix needed for the rgl.points (this is simply a rearrangement to get everything in one table)
  indices <- expand.grid(1:gridRes, 1:gridRes, 1:gridRes)
  den_coor <- expand.grid(dens$x, dens$y, dens$z)
  convertFunc <- function(x, dd){dd[x[1], x[2], x[3]]}
  den_coor$dens <-  apply(indices, 1, convertFunc, dd=dens$d)
  ### Colour palette:
  rbPal <- grDevices::colorRampPalette(c('blue', 'red'))
  ### Define progression of the colour pallette:
  den_coor$col <- rbPal(colres)[as.numeric(cut(den_coor$dens,breaks = colres))]


  #*********************************************************************************
  #   PLOT POINTS/DENSITY   ####
  #*********************************************************************************
  ### Set axis-length if not provided:
  if(is.null(axesLeng)){axesLeng <- c(max(dat$x), max(dat$y), max(dat$z))*1.5}else{axesLeng <- rep(axesLeng, length.out=3)}
  ### Plot empty coordinate system (only if add is FALSE):
  if(!add){
    .coorSystem3D(axesNames, gridCol, axesLeng, axesTicks, zoom)
  }
  ### Start plotting:
  if(grepl('P', whatToPlot[1])){
    rgl::rgl.points(dat[,1:3], color=dat$pointC, size=pointSize, alpha=pointTrans)
  }
  if (grepl('D', whatToPlot[1])){
    rgl::rgl.points(den_coor[,1:3], color=den_coor$col, size=DpointSize, alpha=den_coor$dens/max(den_coor$dens)*Dtransp_fac)   # Here also the progression of the transparency is defined
  }


  #*********************************************************************************
  #   PLOT LINEAR FIT   ####
  #*********************************************************************************
  ### Check what to plot:
  if(plotFit[1]!='no'){
    ### Prepare grid for plotting:
    x.pred <- seq(min(dat$x), max(dat$x), length.out = 100)
    z.pred <- seq(min(dat$z), max(dat$z), length.out = 100)
    xz.grid <- expand.grid("x"=x.pred, "z"=z.pred)
    ### Iterate through the catVar levels:
    for(i in 1:nlevels(dat$catV)){

      #****************************************
      ### In case there is only 1 catVar level:
      if(nlevels(dat$catV)==1){
        if(plotFit[1]=='lin'){
          ### Fit model for coefficients:
          lm.y <- lm(y ~ x + z, dat)
          ### Predict y values:
          y.pred <- matrix(predict(lm.y, newdata = data.frame("x"=xz.grid$x, "z"=xz.grid$z)), ncol = length(x.pred))

        }else if(plotFit[1]=='int'){
          ### Fit model for coefficients:
          lm.y <- lm(y ~ x * z, dat)
          ### Predict y values:
          y.pred <- matrix(predict(lm.y, newdata = data.frame("x"=xz.grid$x, "z"=xz.grid$z)), ncol = length(x.pred))
        }
        #****************************************

      }else if(plotFit[1]=='lin'){
        ### Fit model for coefficients:
        lm.y <- lm(y ~ x + z + catV, dat)
        ### Predict y values:
        y.pred <- matrix(predict(lm.y, newdata = data.frame("x"=xz.grid$x, "z"=xz.grid$z, "catV"=levels(dat$catV)[i])), ncol = length(x.pred))

      }else if(plotFit[1]=='int'){
        ### Fit model for coefficients:
        lm.y <- lm(y ~ x * z + catV, dat)
        ### Predict y values:
        y.pred <- matrix(predict(lm.y, newdata = data.frame("x"=xz.grid$x, "z"=xz.grid$z, "catV"=levels(dat$catV)[i])), ncol = length(x.pred))

      }else if(plotFit[1]=='int2'){
        ### Fit model for coefficients:
        lm.y <- lm(y ~ x + z + catV + catV:x + catV:z, dat)
        ### Predict y values:
        y.pred <- matrix(predict(lm.y, newdata = data.frame("x"=xz.grid$x, "z"=xz.grid$z, "catV"=levels(dat$catV)[i])), ncol = length(x.pred))

      }else if(plotFit[1]=='int3'){
        ### Fit model for coefficients:
        lm.y <- lm(y ~ x * z * catV, dat)
        ### Predict y values:
        y.pred <- matrix(predict(lm.y, newdata = data.frame("x"=xz.grid$x, "z"=xz.grid$z, "catV"=levels(dat$catV)[i])), ncol = length(x.pred))
      }

      ### Set surface color:
      srfCl <- ifelse(test = nlevels(dat$catV) > 1, yes = i, no = "steelblue")
      ### Plot surface of fit:
      rgl::rgl.surface(x.pred, z.pred, y.pred, color = srfCl,
                       alpha = 0.5, lit = FALSE)
    }
  }
}

### Function to plot empty 3d coordinate system (only intended to be used inside nice3DPlot):
.coorSystem3D <- function(axesNames, gridCol, axesLeng, axesTicks, zoom) {
  rgl::rgl.clear()
  rgl::rgl.bg(color = 'white')
  rgl::rgl.lines(c(0, axesLeng[1]), c(0, 0), c(0, 0), color = "black")
  rgl::rgl.lines(c(0, 0), c(0,axesLeng[2]), c(0, 0), color = "black")
  rgl::rgl.lines(c(0, 0), c(0, 0), c(0,axesLeng[3]), color = "black")
  axes <- rbind(c(axesLeng[1], 0, 0), c(0, axesLeng[2], 0),
                c(0, 0, axesLeng[3]))
  rgl::rgl.texts(axes, text = axesNames[1:3], color='black', adj = c(0.5, -0.8))
  if(axesTicks){
    rgl::axis3d(edge = 'x', pos = c(NA,0,0), color='black')
    rgl::axis3d(edge = 'y', pos = c(0,NA,0), color='black')
    rgl::axis3d(edge = 'z', pos = c(0,0,NA), color='black')
  }
  rgl::rgl.viewpoint(theta = -120, phi = 0, zoom=zoom)
  ### grid on floor:
  x <- seq(0,axesLeng[1], length.out = 10)
  z <- seq(0,axesLeng[3],length.out = 10)
  y <- matrix(rep(0, 100), ncol = 10)
  rgl::rgl.surface(x, z, y,
                   color = gridCol,
                   alpha = 0.5, lit = FALSE,front = "lines", back = "lines")
}
