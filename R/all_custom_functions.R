
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
niceUnivPlot <- function(numVar, catVar=NULL, pairedVar=NULL, violin=TRUE, fxdCol=NULL,
                         showMean=TRUE, plot.points=TRUE, bw='nrd0', jitFactor=0.2,
                         add.ylim=0, ylim.cust=NULL, xlab=NULL, ylab=NULL, densScl=0.5,
                         main=NULL, sigGroup=FALSE, sigMu=NULL, multCmp=FALSE,
                         pairCol=NULL, add.lgnd=FALSE, add=FALSE, lnk.means=NULL,
                         lnk.means.lwd=2, pair.lwd=2, ...){


  #*********************************************************************************
  #   CHECK REQUIREMENTS AND PREPARE SOME ARGUMENTS   ####
  #*********************************************************************************
  ### Set important names:
  catVar.nm <- ''
  numVar.nm <- deparse(substitute(numVar))
  ### Check if input is matrix-like object:
  if(inherits(numVar, what = c('matrix', 'data.frame'))){
    ### Check some requirements:
    if(!is.null(catVar)){stop("catVar should not be specified when numVar is supplied as a table")}
    ### Remove non-numeric columns:
    fac_id <- sapply(numVar, inherits, what=c('character', 'factor'))
    if(any(fac_id)){
      numVar <- numVar[, !fac_id, drop=FALSE]
      warning(paste0(sum(fac_id), ' categorical variables were removed for plotting.'))
    }
    ### Store colnames:
    colnms <- colnames(numVar)
    ### Check if pairedVar is supplied:
    if(!is.null(pairedVar)){
      stopifnot(length(pairedVar)==nrow(numVar))   # Must be same length
      ### Get into wide format:
      laps <- lapply(colnms, function(x){data.frame(prVar=pairedVar, value=numVar[,x])})
      ### Extract pairedVar:
      pairedVar <- do.call(rbind, laps)$prVar
      ### Check if pairedCol is supplied as vector:
      if(length(pairCol) > 1){
        ### Check if correct length:
        stopifnot(length(pairCol)==nrow(numVar))
        ### Get into wide format:
        laps <- lapply(colnms, function(x){data.frame(prCol=pairCol, value=numVar[,x])})
        ### Extract pairCol:
        pairCol <- do.call(rbind, laps)$prCol
      }
    }
    ### Get to "wide" format:
    laps <- lapply(colnms, function(x){data.frame(value=numVar[,x], varNm=x)})
    ### Extract num and catVar:
    numVar <- do.call(rbind, laps)$value
    catVar <- factor(do.call(rbind, laps)$varNm, levels=colnms)
  } else{   # In case not a table is supplied as numVar
    ### Check some requirements:
    stopifnot(is.numeric(numVar))
    stopifnot(is.null(ylim.cust) | length(ylim.cust)==2)   # Check valid entries for ylim.cust
    ### Generate catvar:
    catVar.nm <- deparse(substitute(catVar))   # Get name of catVar
    if(is.null(catVar)){
      catVar <- factor(rep(1,length(numVar)))   # 1-level factor
    }else{
      stopifnot(length(catVar)==length(numVar))   # Make sure lengths are same
      catVar <- factor(catVar)   # Dont use as.factor, otherwise non-present levels are not removed.
    }
  }
  ### Set ylim:
  ylms <- c(min(numVar, na.rm = TRUE), max(numVar, na.rm = TRUE))
  ylms <- ylms + c(-add.ylim, add.ylim)   # Add the add.ylim


  #*********************************************************************************
  #   RUN GROUP COMPARISONS   ####
  #*********************************************************************************
  ### Only run if catVar has entries of multiple levels:
  drawGrCmp <- FALSE   # Indicator if group comparisons can be plotted
  if(nlevels(catVar) > 1 & sigGroup){
    ### List group comparisons in table:
    grInd <- data.frame(t(utils::combn(1:nlevels(catVar),2)))
    ### Apply wilcoxon tests to compare corresponding groups (collect pvalues):
    if(is.null(pairedVar)){   # Unpaired test
      wilxP <- apply(grInd, 1, function(x){ wilcox.test(x = numVar[as.numeric(catVar)==x[1]],
                                                        y = numVar[as.numeric(catVar)==x[2]])$p.value })
    }else{   # Paired test
      ### Check whether length of pairedVar is correct:
      stopifnot(length(pairedVar)==length(numVar), is.factor(pairedVar))
      ### Function to perform correct pairwise test (have to be careful that same cases are paired)
      pairSel <- function(x){
        dd <- data.frame(numVar, catVar, pairedVar)
        xx <- dd[as.numeric(dd$catVar)==x[1],]
        yy <- dd[as.numeric(dd$catVar)==x[2],]
        ### Check which cases are present in both groups:
        lvs <- 1:nlevels(pairedVar)
        lvsOk <- which(lvs %in% as.numeric(xx$pairedVar) & lvs %in% as.numeric(yy$pairedVar))
        xx.1 <- xx[as.numeric(xx$pairedVar) %in% lvsOk,]
        yy.1 <- yy[as.numeric(yy$pairedVar) %in% lvsOk,]
        ### Apply wilcoxon test:
        wilcox.test(x = xx.1[order(xx.1$pairedVar),'numVar'],
                    y = yy.1[order(yy.1$pairedVar),'numVar'], paired = TRUE)$p.value
      }
      wilxP <- apply(grInd, 1, pairSel)
    }
    ### Combine to table (and apply Bonferonni correction):
    grCmp.0 <- cbind(grInd, wilxP)
    if(multCmp){grCmp.0$wilxP <- wilxP*nrow(grInd)}   # Apply Bonferroni correction
    ### Add stars to table:
    grCmp.0$strs <- symnum(grCmp.0$wilxP, corr = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, ifelse(max(grCmp.0$wilxP) > 1, max(grCmp.0$wilxP)+1, 1)),
                           symbols = c("***", "**", "*", ""), legend = FALSE)
    ### Filter out only significant tests:
    grCmp <- grCmp.0[grCmp.0$wilxP < 0.05,]
    ### Only proceed if there are significant tests:
    if(nrow(grCmp) > 0){
      ### Add y-values for plot (should not be too far apart):
      grCmp$yVal <- max(ylms) + (abs(min(ylms)-max(ylms))*0.05 * 1:nrow(grCmp))
      ### Adapt ylms:
      ylms <- c(ylms[1], max(grCmp$yVal))
      ### Set indicator to plot to TRUE:
      drawGrCmp <- TRUE
    }else{warning('There were no significant group differences.')}
  }


  #*********************************************************************************
  #   RUN MU COMPARISONS   ####
  #*********************************************************************************
  ### Test each group against a mean value Mu using a wilcoxon test:
  drawMu <- FALSE   # Indicator to draw results in plot
  if(!is.null(sigMu)){
    ### Iterate through groups:
    p.mu <- NA
    for(i in 1:nlevels(catVar)){
      dmu.i <- numVar[as.numeric(catVar)==i]
      p.mu[i] <- wilcox.test(dmu.i, mu = sigMu)$p.value
    }
    # Apply Bonferroni correction:
    if(multCmp){p.mu <- p.mu*nlevels(catVar)}
    ### Create table:
    d.mu <- data.frame("categ"=1:nlevels(catVar), "pmu"=p.mu)
    ### Add stars:
    d.mu$strs <- symnum(d.mu$pmu, corr = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, ifelse(max(d.mu$pmu) > 1, max(d.mu$pmu)+1, 1)),
                        symbols = c("***", "**", "*", ""), legend = FALSE)
    ### Filter out only significant tests:
    d.mu <- d.mu[d.mu$pmu < 0.05,]
    ### Only proceed if there are significant tests:
    if(nrow(d.mu) > 0){
      ### Adapt ylms:
      if(sigMu < ylms[1]){
        ylms[1] <- sigMu
      }else if(sigMu > ylms[2]){
        ylms[2] <- sigMu
      }
      ### Set indicator to plot to TRUE:
      drawMu <- TRUE
    }else{warning('There were no significant differences from tested location Mu.')}
  }


  #*********************************************************************************
  #   PLOT POINTS   ####
  #*********************************************************************************
  ### Get title name:
  main.nm <- deparse(substitute(main))
  ### Plugin custom y-values:
  if(!is.null(ylim.cust)){
    ylms <- ylim.cust
  }
  ### Plugin custom axes labels:
  if(is.null(xlab)){
    xlab <- ifelse(catVar.nm=='NULL', '', catVar.nm)
  }
  if(is.null(ylab)){
    ylab <- numVar.nm
  }
  ### Start plot (empty):
  if(!add){
    plot(x=1,y=1,
         type='n',
         xlim = c(0,nlevels(catVar)+1),
         ylim = ylms,
         xaxt = 'n',
         ylab = ylab,
         xlab = xlab,
         main = ifelse(main.nm=='NULL', paste0(numVar.nm, ' Plot'), main))
  }
  ### Add points:
  if(plot.points){
    points(x = jitter(as.numeric(catVar), factor = jitFactor),
           y = numVar,
           col = if(is.null(fxdCol)){catVar}else{fxdCol},
           ...)
  }
  ### Add x-axis labels:
  if(!catVar.nm=='NULL' & !add){
    axis(1, at = 1:nlevels(catVar), labels = levels(catVar))
  }
  ### Add legend:
  if(nlevels(catVar) > 1 & add.lgnd){   # Only needed with multiple levels
    legend('bottomright', legend = levels(catVar),
           pch=1,
           col=if(is.null(fxdCol)){1:nlevels(catVar)}else{fxdCol})
  }


  #*********************************************************************************
  #   PLOT LINES IN PAIRED CASE   ####
  #*********************************************************************************
  ### Only run if a "pairedVar" factor was provided:
  if(!is.null(pairedVar) & nlevels(catVar) > 1){
    ### Make sure there are no levels with no entries:
    pairedVar <- factor(pairedVar)
    ### Check if every case has maximally one entry per catVar level:
    if(max(table(catVar, pairedVar)) > 1){ warning('There are cases with multiple entries for a catVar level.') }
    ### Define colour of lines:
    if(is.null(pairCol)){
      pCol <- as.numeric(pairedVar)   # Each case gets its own colour
    } else if(length(pairCol)==length(pairedVar)){
      if(is.factor(pairCol)){
        pCol <- as.numeric(pairCol)   # Colour of lines according to pairCol factor
      }
      if(is.numeric(pairCol)){
        br_ramp <- colorRampPalette(c('red','blue'))
        pCol <- br_ramp(10)[as.numeric(cut(pairCol, breaks=10))]   # Colour ranging from blue to red according to value of pairCol numeric
      }
    }else{pCol <- pairCol}   # When pairCol is one fixed colour
    ### Draw lines per case in for-loop:
    dd.0 <- data.frame(numVar, catVar, pairedVar, pCol)
    for(i in 1:nlevels(pairedVar)){
      dd <- dd.0[as.numeric(dd.0$pairedVar)==i,]   # Select entries of case i
      dd <- dd[order(dd$catVar), ]   # Make sure order of catVar is correct
      ### Check whether pCol values are all equal for case:
      if(length(unique(dd$pCol)) != 1){warning('There are cases with different values of the pairCol factor.')}
      ### Plot line:
      lines(x = dd$catVar, y = dd$numVar, col=dd$pCol[1], lwd=pair.lwd)
    }
  }


  #*********************************************************************************
  #   ADD VIOLIN LINES   ####
  #*********************************************************************************
  ### Add the violin lines:
  if(violin){
    L <- list()
    for(i in 1:nlevels(catVar)){
      d <- density(numVar[as.numeric(catVar)==i], na.rm = TRUE, bw = bw)
      ### Get min and max value of numVar:
      minNum.i <- min(numVar[as.numeric(catVar)==i], na.rm = TRUE)
      maxNum.i <- max(numVar[as.numeric(catVar)==i], na.rm = TRUE)
      ### Remove density values falling outside the min-max range:
      denx <- d$x[d$x > minNum.i & d$x < maxNum.i]
      deny <- d$y[d$x > minNum.i & d$x < maxNum.i]
      ### Turn density values at each end to zero to cut off the density curve:
      deny[1] <- 0
      deny[length(deny)] <- 0
      L[[i]] <- data.frame(xd=denx, yd=deny)
    }
    names(L) <- levels(catVar)
    ### We have to scale the densities, need the maximum value for that:
    maxD <- max(do.call(c, lapply(L, function(x){x$yd})))
    cexD <- densScl/maxD
    ### Now plot the densities:
    for(i in 1:nlevels(catVar)){
      lines(L[[i]]$yd*cexD + i, L[[i]]$xd, col= if(is.null(fxdCol)){i}else{fxdCol}, lwd=3)
      lines((-L[[i]]$yd)*cexD + i, L[[i]]$xd, col=if(is.null(fxdCol)){i}else{fxdCol}, lwd=3)
    }
  }


  #*********************************************************************************
  #   ADD MEAN LINES AND CONNECTIONS   ####
  #*********************************************************************************
  ### Add the mean-value lines:
  if(showMean){
    for(i in 1:nlevels(catVar)){
      mVal <-  mean(numVar[as.numeric(catVar)==i], na.rm = TRUE)
      segments(x0 = i-0.3, y0 = mVal, x1 = i+0.3, y1 = mVal, col = if(is.null(fxdCol)){i}else{fxdCol}, lwd = 3)
    }
  }
  ### Add mean connections:
  if(!is.null(lnk.means)){
    mVal <-  tapply(numVar, INDEX = catVar, FUN = mean, na.rm=TRUE)
    lines(x = 1:nlevels(catVar), y = mVal, col=lnk.means, lwd=lnk.means.lwd)
  }


  #*********************************************************************************
  #   ADD LINES OF GROUP COMPARISONS   ####
  #*********************************************************************************
  ### Only run if there are significant results:
  if(drawGrCmp){
    ### Draw lines:
    for(i in 1:nrow(grCmp)){
      xx <- as.numeric(grCmp[i,1:2])
      yy <- rep(grCmp[i,'yVal'], 2)
      lines(xx, yy)
      ### Find a good ticklength:
      tickL <- (abs(min(ylms)-max(ylms))*0.01)
      lines(c(xx[1], xx[1]), c(yy[1], yy[1] - tickL))
      lines(c(xx[2], xx[2]), c(yy[2], yy[2] - tickL))
      ### Add stars:
      text(x = mean(xx), y = yy[1]+tickL, labels = grCmp[i,'strs'])
    }
  }


  #*********************************************************************************
  #   ADD LINES OF MU COMPARISONS   ####
  #*********************************************************************************
  ### Only run if there are significant results:
  if(drawMu){
    ### Draw lines:
    for(i in 1:nrow(d.mu)){
      xx <- rep(d.mu$categ[i]-0.4, 2)
      yy <- c(mean(numVar[as.numeric(catVar)==d.mu$categ[i]], na.rm = TRUE), sigMu)
      lines(xx,yy)
      ### Find good ticklength:
      tickL <- nlevels(catVar)/100
      lines(c(xx[1], xx[1] + tickL), c(yy[1], yy[1]))
      lines(c(xx[2], xx[2] + tickL), c(yy[2], yy[2]))
      ### Add stars:
      text(x = xx[1]-2*tickL, y = mean(yy), labels = d.mu[i,'strs'])
      ### Draw Mu value as line:
      abline(h = sigMu, lty=2)
    }
  }
}


#*********************************************************************************
#   NICE PAIRS PLOT   ####
#*********************************************************************************
nicePairsPlot <- function(x, catVar = NULL, breaks = "Sturges", density = FALSE, jitter = FALSE, jitFactor = 1, loess = FALSE, swtchPan = FALSE, exclude = c("none", "numeric", "factor"), keepOrder = FALSE, facsAtBegin = FALSE, cex.diag = 2, cex.offdiag = 1, cex.mean = 2, cex.lvls = 1){


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   PREPARE DATA   ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Make sure x is a data frame:
  x <- data.frame(x)
  ### Check some things about exclude option:
  if(any(exclude %in% colnames(x))){   # In case of manually selecting vars
    if(!all(exclude %in% colnames(x))){warning('Not all variables to exclude are present in the data')}
    if(any(exclude %in% as.character(formals(nicePairsPlot)$exclude)[-1])) warning(paste0('The exclude argument is appearing in the colnames but also includes one of the general exclusion options (',  paste(as.character(formals(nicePairsPlot)$exclude)[-1], collapse = ', '), ')'))
  }else{
    exclude <- match.arg(exclude)
    cls_toex <- switch (exclude,
                        none = '',
                        numeric = c('numeric', 'integer'),
                        factor = c('character', 'factor')
    )
    ### Get var names to exclude:
    ex_log <- sapply(x, inherits, what=cls_toex)
    exclude <- colnames(x)[ex_log]
  }
  ### Check some things regarding catVar:
  if(!is.null(catVar)){
    stopifnot(nrow(x)==length(catVar))
    stopifnot(inherits(catVar, 'factor'))
  }
  ### Remove variables to exclude:
  x <- x[, !(colnames(x)%in%exclude), drop=FALSE]
  ### Turn every character to factor:
  lfac <- lapply(x, function(a){
    if(inherits(a, 'character')) {rval <- factor(a)}else{rval <- a}
  })
  x <- do.call(data.frame, lfac)
  ### Reorder columns to have factors together:
  if(!keepOrder){
    fac_id <- sapply(x, inherits, what='factor')   # Indices
    if(!facsAtBegin){
      x <- cbind(x[, !fac_id, drop=FALSE], x[, fac_id, drop=FALSE])   # At end
    }else{
      x <- cbind(x[, fac_id, drop=FALSE], x[, !fac_id, drop=FALSE])   # At beginning
    }
  }
  ### Store colnames:
  nms <- colnames(x)
  ### Store whether factor:
  fac_id <- sapply(x, inherits, what='factor')
  ### Make sure data has multiple variables:
  if(ncol(x) < 2){stop('There must be at least two variables for plotting.')}


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   DIFFERENT PLOT FUNCTIONS   ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Barplot factors:
  barpl_fac <- function(vx, vy){
    ### Only one variable needed:
    v <- vx
    ### xylims:
    xlims <- c(0, nlevels(v)) + 0.5
    ylims <- c(0, 1.4*max(table(v)))
    ### Prepare Plot:
    plot(x = NA, xlab='', ylab='', yaxt='n', xaxt='n', xlim=xlims, ylim=ylims, yaxs='i', xaxs='i')
    ### Background:
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#b6f0fc")
    ### Add rectangles:
    bw <- 0.43   # bar width
    recpos <- 1:nlevels(v)
    rect(recpos - bw, 0, recpos + bw, table(v), col = 'cyan')
    ### Axes:
    if(edge_rw!=0) axis(edge_cl)
    ### Level names:
    if(nlevels(v) < 10){
      text(x = recpos, y = ylims[1]+0.05*ylims[2], labels = levels(v), srt=90, pos=4, offset=0, cex=cex.lvls)
    }
    ### Title:
    text(x = mean(par()$usr[1:2]), y = ylims[2]*0.85, labels = nms[rw], cex=cex.diag)
  }
  #************************************
  ### Histogramm:
  histgr <- function(vx, vy){
    ### Only one variable needed:
    v <- vx
    ### Histogramm data:
    h <- hist(v, breaks = breaks, plot = FALSE)
    hbrks <- h$breaks
    nB <- length(hbrks)
    y <- h$counts
    ### xylims:
    xlims <- range(hbrks)
    ylims <- c(0, 1.4*max(y))
    ### Prepare plot:
    plot(x = NA, xlab='', ylab='', yaxt='n', xaxt='n', xlim=xlims, ylim=ylims, yaxs='i', xaxs='i')
    ### Add axes:
    if(edge_cl!=0) axis(edge_cl)
    ### Add rectangles:
    rect(hbrks[-nB], 0, hbrks[-1], y, col = 'cyan')
    ### Title:
    text(x = mean(par()$usr[1:2]), y = ylims[2]*0.85, labels = nms[rw], cex=cex.diag)
    ### Add density:
    if(density){
      tryd <- try(d <- density(vx, na.rm = TRUE), silent = TRUE)
      if (class(tryd) != 'try-error') {
        d$y <- d$y/max(d$y)*max(y)
        lines(d, col='darkblue')
      }
    }
  }
  #************************************
  ### Scatterplot:
  scattpl <- function(vx, vy){
    ### Indicator if mean vals should be plotted:
    pltm <- FALSE
    ### xylims and jitter:
    jitx <- jity <- 0   # Default
    jitf_fct <- 0.5   # Have jitter a bit smaller for factors
    ### Var x:
    if(fac_id[cl]){
      xlims <- c(0, nlevels(vx)) + 0.5
      jitx <- jitter(as.numeric(vx), factor=jitFactor*jitf_fct) - as.numeric(vx)   # Create jitter
    }else{
      xlims <- range(vx, na.rm = TRUE) + (diff(range(vx, na.rm = TRUE))*0.05)*c(-1,1)
      jitx <- if (jitter){ jitter(as.numeric(vx), factor=jitFactor) - as.numeric(vx)}else{0}
    }
    ### Var y:
    if(fac_id[rw]){
      ylims <- c(0, nlevels(vy)) + 0.5
      jity <- jitter(as.numeric(vy), factor=jitFactor*jitf_fct) - as.numeric(vy)   # Create jitter
    }else{
      ylims <- range(vy, na.rm = TRUE) + (diff(range(vy, na.rm = TRUE))*0.05)*c(-1,1)
      jity <- if (jitter){ jitter(as.numeric(vy), factor=jitFactor) - as.numeric(vy)}else{0}
    }
    ### Prepare plot:
    plot(x = NA, xlab='', ylab='', yaxt='n', xaxt='n', xlim=xlims, ylim=ylims, yaxs='i', xaxs='i')
    ### Background:
    if(all(fac_id[c(rw, cl)])){
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#b6f0fc")
    }else if(any(fac_id[c(rw, cl)])){
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#e8fbff")
      ### Indicator if mean vals should be plotted:
      pltm <- TRUE
    }
    ### Add points:
    cls <- if(is.null(catVar)){1}else{catVar}   # Set colour
    points(x = as.numeric(vx) + jitx, y = as.numeric(vy) + jity, col=cls)
    ### Add mean values (only for plots with numeric and factor):
    if(pltm){
      ### Check which one is factor:
      if(fac_id[rw]){
        points(x = tapply(vx, vy, mean, na.rm=TRUE), y = 1:nlevels(vy), col='cyan', cex=cex.mean, pch=19)
      }else{
        points(x = 1:nlevels(vx), y = tapply(vy, vx, mean, na.rm=TRUE), col='cyan', cex=cex.mean, pch=19)
      }
    }
    ### Add axes:
    if(edge_rw!=0){
      if(fac_id[cl]){
        axis(edge_rw, at = 1:nlevels(vx), labels = levels(vx))
      }else{
        axis(edge_rw)
      }
    }
    if(edge_cl!=0){
      if(fac_id[rw]){
        axis(edge_cl, at = 1:nlevels(vy), labels = levels(vy))
      }else{
        axis(edge_cl)
      }
    }
    ### Loess line (only for numerics):
    if(loess & !fac_id[rw] & !fac_id[cl]){
      tryd <- try(lml <- suppressWarnings(loess(vy ~ vx, degree = 1, family = "symmetric")), silent=TRUE)
      if(class(tryd)!='try-error'){   # In case of no error
        tempx <- data.frame(vx = seq(min(vx, na.rm = TRUE),
                                     max(vx, na.rm = TRUE), length.out = 50))
        pred <- predict(lml, newdata = tempx)
        lines(x=tempx$vx, y=pred, col='cyan', lty=2, lwd=2)
      }
    }
  }


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   DIFFERENT STATISTICAL TEST FUNCTIONS   ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Simple correlation:
  cortst <- function(vx, vy){
    ### Add try statement because of potential errors:
    tryd <- try({
      cort <- cor.test(vx, vy)
      p <- cort$p.value
      psymbs <- c("***", "**", "*", "")
      star <- symnum(p, corr = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = psymbs, legend = FALSE)
      txt <- paste0('r=   ', round(cort$estimate, digits=2), star)
      ### Size of numbers:
      poscex <- seq(from=1, by=0.8, length.out=4)   # Possible sizes
      poscex <- poscex[length(poscex):1]   # Reverse
      cex <- poscex[psymbs %in% star] * cex.offdiag
    }, silent = TRUE)
    ### In case of error:
    if(class(tryd)=='try-error'){
      txt <- 'Err'
      cex <- 1
    }
    ### Write statistic and add stars:
    plot(x = NA, xlab='', ylab='', yaxt='n', xaxt='n', xlim=0:1, ylim=0:1, yaxs='i', xaxs='i')
    text(0.5, 0.5, txt, cex = cex)
  }
  #************************************
  ### ANOVA:
  aovtst <- function(vx, vy){
    ### Which var is factor:
    if(fac_id[rw]){ yaov <- vx; xaov <- vy }else{ yaov <- vy; xaov <- vx }
    ### Add try statement because of potential errors:
    tryd <- try({
      res <- anova(lm(yaov ~ xaov))
      fval <- round(res$`F value`[1], 1)
      p <- res$`Pr(>F)`[1]
      psymbs <- c("***", "**", "*", "")
      star <- symnum(p, corr = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = psymbs, legend = FALSE)
      txt <- paste0('F=   ', fval, star)
      ### Size of numbers:
      poscex <- seq(from=1, by=0.8, length.out=4)   # Possible sizes
      poscex <- poscex[length(poscex):1]   # Reverse
      cex <- poscex[psymbs %in% star] * cex.offdiag
    }, silent = TRUE)
    ### In case of error:
    if(class(tryd)=='try-error'){
      txt <- 'Err'
      cex <- 1
    }
    ### Setup plot:
    plot(x = NA, xlab='', ylab='', yaxt='n', xaxt='n', xlim=0:1, ylim=0:1, yaxs='i', xaxs='i')
    ### Background:
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#e8fbff")
    ### Write statistic and add stars:
    text(0.5, 0.5, txt, cex = cex)
  }
  #************************************
  ### Chisq test:
  chitst <- function(vx, vy){
    ### Add try statement because of potential errors:
    tryd <- try({
      res <- chisq.test(table(vx, vy))
      xval <- round(res$statistic, 1)
      p <- res$p.value
      psymbs <- c("***", "**", "*", "")
      star <- symnum(p, corr = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     symbols = psymbs, legend = FALSE)
      txt <- paste0('Chi=   ', xval, star)
      ### Size of numbers:
      poscex <- seq(from=1, by=0.8, length.out=4)   # Possible sizes
      poscex <- poscex[length(poscex):1]   # Reverse
      cex <- poscex[psymbs %in% star] * cex.offdiag
    }, silent = TRUE)
    ### In case of error:
    if(class(tryd)=='try-error'){
      txt <- 'Err'
      cex <- 1
    }
    ### Prepare plot:
    plot(x = NA, xlab='', ylab='', yaxt='n', xaxt='n', xlim=0:1, ylim=0:1, yaxs='i', xaxs='i')
    ### Background:
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#b6f0fc")
    ### Write statistic and add stars:
    text(0.5, 0.5, txt, cex = cex)
  }


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   CREATE PLOT INDICATION MATRIX   ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Visual plots:
  ### Keys:
  pltnms <- c('histgr', 'barpl_fac', 'scattpl')
  ### Matrix:
  plt_m <- matrix(3, nrow=ncol(x), ncol = ncol(x))
  ### Diagonal:
  diag(plt_m) <- fac_id+1
  ### Turn to text matrix:
  plt_m_chr <- matrix(pltnms[as.vector(plt_m)], nrow=ncol(x))
  ### Test matrix:
  ### Keys:
  tstnms <- c('cortst', 'aovtst', 'chitst')
  ### Matrix:
  tst_m <- outer(fac_id, fac_id, "+")+1
  ### Turn to text matrix:
  tst_m_chr <- matrix(tstnms[as.vector(tst_m)], nrow=ncol(x))
  ### Combine the two matrices:
  if(!swtchPan){
    plt_m_chr[lower.tri(plt_m_chr)] <-  tst_m_chr[lower.tri(tst_m_chr)]
  }else{
    plt_m_chr[upper.tri(plt_m_chr)] <-  tst_m_chr[upper.tri(tst_m_chr)]
  }


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #   CREATE PLOT   ####
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ### Prepare plot:
  olpar <- parSave()
  on.exit(parReset(olpar))   # Make sure par is reset
  par(mfrow=c(ncol(x), ncol(x)))
  ### Margins:
  mars <- rep(0.5, 4)
  par(mar=mars)
  par(oma = rep(3, 4))
  ### Loop through plot grid:
  for(rw in 1:ncol(x)){
    ### Check if on edge:
    edge_rw <- if(rw==1){3} else if(rw==ncol(x)){1} else{0}
    for(cl in 1:ncol(x)){
      ### Check if on edge:
      edge_cl <- if(cl==1){2} else if(cl==ncol(x)){4} else{0}
      ### Apply correct plot:
      do.call(plt_m_chr[rw, cl], list(vx = x[,cl], vy = x[,rw]))
    }
  }
}

#*********************************************************************************
#   NICE 3D PLOT   ####
#*********************************************************************************
nice3DPlot <- function(X = NULL, whatToPlot = c('P','D','PD'), plotFit = c('no','lin','int','int2','int3'), catVar = factor(1), covMat = NULL, means = NULL, nSim = 500, pointCol = 1, colRamp = NULL, pointSize = NULL, spheres = FALSE, pointTrans = 0.8, Dtransp_fac = 0.06, colres = 50, gridRes = 30, h = 0.3, DpointSize = 30, axesNames = NULL, axesTicks = FALSE, gridCol = 'grey', axesLeng = NULL, zoom = 1, add = FALSE, htmlFilename=NULL, ...){
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
    ### In case colour ramp up is asked for:
    if(!is.null(colRamp)){
      ### Some checks:
      if(any(is.na(colRamp))){warning('There are missings in the colRamp variable')}
      if(length(colRamp) != nrow(X)){warning('colRamp does not have the same length as the data to be plotted.')}
      ### Prepare palette:
      ncl_pl <- 100
      col_pl <- colorRampPalette(c('navy', 'hotpink'))(ncl_pl)
      ### Cut likel. values into ranges:
      col_indx <- as.numeric(cut(colRamp, breaks=ncl_pl))
      ### Point colour vector:
      pointCol <- col_pl[col_indx]
    }
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
    #*********************
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
    if(spheres){   # Should spheres be plotted or points
      if(is.null(pointSize)){pointSize <- 0.05}
      rgl::rgl.spheres(dat[,1:3], color=dat$pointC, radius=pointSize, alpha=pointTrans)
    }else{
      if(is.null(pointSize)){pointSize <- 5}
      rgl::rgl.points(dat[,1:3], color=dat$pointC, size=pointSize, alpha=pointTrans, point_antialias = TRUE)   # point_antialias makes points round in WebGL (or use rgl.spheres alternatively)
    }
  }
  if (grepl('D', whatToPlot[1])){   # Plot density
    rgl::rgl.points(den_coor[,1:3], color=den_coor$col, size=DpointSize, point_antialias = TRUE,
                    alpha=den_coor$dens/max(den_coor$dens)*Dtransp_fac)   # Here also the progression of the transparency is defined
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

  #*********************************************************************************
  #   SAVE HTML FILE   ####
  #*********************************************************************************
  ### Check if filename exists:
  if(!is.null(htmlFilename)){
    ### Create widget:
    widget <- rgl::rglwidget(...)
    ### Save as html:
    htmlwidgets::saveWidget(widget, htmlFilename)
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


#*********************************************************************************
#   MIXED MODEL DGP   ####
#*********************************************************************************
mmdgp <- function(n=200, nC=5, sd_S=3, sd_C=5, sd_e=2, b0=50, tb=c(0,4,10,4,0),
                  genb=10, ageb=1){
  ### Check some condition:
  if(n%%nC!=0){stop('There has to be an equal number of people in each class. Currently, the number of people is not divisible by the number of classes.')}
  ### Generate subjects:
  id <- factor(paste0('S', 1:n), levels = paste0('S', 1:n))
  ### Generate random intercept of subjects:
  id.ri <- rnorm(n, sd=sd_S)
  ### Generate age variable:
  age <- ceiling(runif(n, min = 20, max = 70))
  ### Generate class factor:
  class <- factor(paste0('Class', 1:nC), levels = paste0('Class', 1:nC))
  ### Generate class random intercept:
  clss.ri <- rnorm(nC, sd=sd_C)
  ### Generate gender variable:
  gen <- 0:1
  ### Combine to data frame:
  ds <- data.frame(id, class, gen, age, id.ri, clss.ri)
  ds <- ds[order(ds$class, ds$id),]   # Order after classes and individuals
  ds$gen <- sample(ds$gen)   # mix up gender
  ### Generate target variable for each day:
  ### First replicate each row in data frame for each day:
  dat.0 <- replicate(length(tb), ds, simplify = FALSE)
  dat.0 <- Reduce(f = function(a,b)rbind(a,b), x = dat.0)
  dat.0 <- dat.0[order(dat.0$class, dat.0$id),]   # Reorder
  ### Add day variable:
  dat.0$day <- 1:length(tb)
  ### Add random error:
  dat.0$err <- rnorm(nrow(dat.0), sd=sd_e)
  ### Generate target variable for each day (bit messy...):
  dd <- dat.0[, -which(colnames(dat.0) %in% c('id', 'class', 'day'))]   # Remove some variables
  dd$dayEff <- tb
  dd$b0 <- b0
  dat.0$y <- apply(dd, 1, function(x){
    t(as.matrix(x))%*%c(genb, ageb, 1, 1, 1, 1, 1)
  })
  ### Remove the columns not needed:
  dat <- dat.0
  dat$id.ri <- NULL
  dat$clss.ri <- NULL
  dat$err <- NULL
  ### Gender to factor:
  dat$gen <- factor(dat$gen, levels = c(0,1), labels = c('M','F'))
  ### Day to factor:
  dat$day <- as.factor(dat$day)
  ### Return object:
  return(dat)
}


#*********************************************************************************
#   WIDE TO LONG DATA FORMAT    ####
#*********************************************************************************
wideToLong <- function(x, nRep=NULL, ind='T*_', indCust=NULL, repColnm='repIdentifier'){
  ### Check some conditions:
  if(is.null(nRep) & is.null(indCust)){
    stop('Either nRep or indCust must be supplied')
  }
  ### Get names of rep-identifier:
  if(is.null(indCust)){
    rnms <- sapply(1:nRep, function(y)gsub(pattern = '\\*', replacement = y, ind))
  }else{
    rnms <- indCust
  }
  ### Check whether the columns are all named correctly:
  tvr <- grep(paste0(rnms, collapse = '|'), colnames(x))   # ColumnNr which are repeated measures
  tvars <- colnames(x[,tvr])   # colnames
  ### Remove the rep-identifyer:
  tvars.v <- sub(pattern = paste0(rnms, collapse = '|'), replacement = '', perl = TRUE, tvars)
  ### Remove the var names:
  tvars.t <- gsub(pattern = paste0('(', paste0(rnms, collapse = '|'),')(*SKIP)(*FAIL)|.'),
                  replacement = '', perl = TRUE, tvars)
  #' This works like: What_I_want_to_avoid(*SKIP)(*FAIL)|What_I_want_to_match
  #' See also these links (should really learn regex):
  #' https://stackoverflow.com/questions/24534782/how-do-skip-or-f-work-on-regex
  #' https://stackoverflow.com/questions/38712946/gsub-everything-except-specified-characters
  ### Data frame to test names:
  tvars.d <- data.frame(tvars.v, tvars.t)
  if(!all(table(tvars.d)==1)){
    print(table(tvars.d))
    print('There is a problem with the naming of the repeated variables. Check returned table (printed above), should all be equal to 1.')
    return(table(tvars.d))
  }
  ### Iterate through every single row and collect long format in list:
  L <- list()
  tvars.vu <- unique(tvars.v)
  for(i in 1:nrow(x)){
    tmp.f <- x[i, -tvr, drop=FALSE]   # Fixed variables
    tmp.v <- x[i, tvr, drop=FALSE]   # Repeated measures variables
    tmp <- tmp.f[rep(1, length(rnms)),, drop=FALSE]   # Repeat row for each repeated measurement
    tmp[, repColnm] <- rnms   # Add rep-identifiers
    for(j in 1:length(tvars.vu)){   # Iterate through the individual tests
      p <- paste0(rnms, tvars.vu[j])
      p.ind <- match(p, colnames(tmp.v))
      tmp[, tvars.vu[j]] <- as.vector(t(tmp.v[1, p.ind]))
    }
    L[[i]] <- tmp
  }
  ### Put everything together:
  d <- do.call(rbind, L)
  ### Turn rep-identifier to factor:
  d[, repColnm] <- factor(d[, repColnm], levels = rnms, labels = rnms)
  ### Return object:
  return(d)
}


#*********************************************************************************
#   NICE NA PLOT    ####
#*********************************************************************************
niceNaPlot <- function(x, IDvar=NULL, show_xlab=TRUE){
  ### If supplied, add IDvar as rownames and then remove it:
  if(is.null(IDvar)){
    IDvar <- paste0('Obs', 1:nrow(x))
    rownames(x) <- IDvar
  }else{
    ### In case there are repetitions of identifier values (rownames must be unique):
    if(length(unique(x[,IDvar])) != nrow(x)){
      x <- x[order(x[,IDvar]),]   # Make sure identifier are ordered
      id.rle <- rle(as.character(x[,IDvar]))
      id.uni <- paste0(rep(id.rle$values, times = id.rle$lengths), "_",
                       unlist(lapply(id.rle$lengths, seq_len)))
    }else{id.uni <- x[,IDvar]}
    ### Add as rownames:
    rownames(x) <- id.uni
    x[,IDvar] <- NULL
  }
  ### Create NA table:
  x <- is.na(x)*1
  ### Apply hierarchical clustering to find good ordering:
  dis <- dist(x)
  hc <- hclust(dis, method = 'ward.D2')
  x <- x[hc$order,]
  ### Create plot:
  olmar <- newmar <- par('mar')
  newmar[c(2,4)] <- newmar[c(2,4)]+4
  par('mar'=newmar)
  xx <- 1:nrow(x)
  yy <- 1:ncol(x)
  image(xx,yy,x,col=c('lightblue', 'darkgrey'), xaxt='n',
        yaxt='n', xlab='', ylab = '')
  ### Set xlab:
  if(show_xlab){
    xlbs <- rownames(x)
  }else{
    xlbs <- rep('', nrow(x))
  }
  ### Add axes:
  axis(side = 1, at = 1:nrow(x), labels = xlbs, las=2)
  axis(side = 2, at = 1:ncol(x), labels = colnames(x), las=1)
  par(xpd=TRUE)
  legend(x = par('usr')[2], y=par('usr')[4], legend = c('present', 'missing'),
         pch=15, col=c('lightblue', 'darkgrey'), pt.cex=2, bg='grey90', box.lty = 'blank')
  par('mar'=olmar)
  par(xpd=FALSE)   # Return to default
  ### In case one wants the ordered is.na-table:
  silentReturn <- x
}

#*********************************************************************************
#   LONG TO WIDE DATA FORMAT    ####
#*********************************************************************************
longToWide <- function(x, IDvar, repColnm, repVars, colNmStr=''){
  ### Make sure IDvar and repColnm are factors:
  stopifnot(is.factor(x[,repColnm]))
  stopifnot(is.factor(x[,IDvar]))
  ### Remove cases with missing values in repColnm or IDvar:
  nasm.rep <- sum(is.na(x[,repColnm]))
  nasm.id <- sum(is.na(x[,IDvar]))
  if(nasm.rep != 0){
    warning(paste0('There are ', nasm.rep,' missing value(s) in the ', repColnm ,' variable. Will remove these observations.'))
    x <- x[!is.na(x[,repColnm]),]
  }
  if(nasm.id != 0){
    warning(paste0('There are ', nasm.id,' missing value(s) in the IDvar variable. Will remove these observations.'))
    x <- x[!is.na(x[,IDvar]),]
  }
  ### Check whether there are multiple entries for a repColnm level for one ID:
  if(max(table(x[,IDvar], x[,repColnm])) > 1){
    print(table(x[,IDvar], x[,repColnm]))
    stop(paste0('There are multiple entries for a ', repColnm, ' level for a level of ', IDvar, '. Check the printed table to find out more.'))
  }
  ### Get the levels of repColnm:
  lvs <- levels(x[,repColnm])
  ### Split data according to observation-ID:
  x.sp <- split(x, x[, IDvar])
  ### Loop through the IDs:
  xwL <- list()
  for(i in 1:length(x.sp)){
    ### Get data:
    d <- x.sp[[i]]
    ### Check whether fixed variables have varying values:
    dfi <- d[, !(colnames(d) %in% repVars) & colnames(d)!=repColnm]
    if(!all(sapply(dfi, function(x){length(unique(x))==1}))){
      warning(paste0('There are varying values in presumably fixed variable(s) for observation-ID ', levels(x[,IDvar])[i],'. Will use the values of the first row.'))
    }
    ### Check whether levels of repColnm are missing and adjust:
    if(length(unique(d[,repColnm])) < length(lvs)){
      ### Add number of missing rows:
      mss <- length(lvs) - length(unique(d[,repColnm]))
      dd <- d[rep(1, mss),]
      ### Fill in appropriate values:
      dd[, repVars] <- NA
      dd[, repColnm] <- lvs[!lvs %in% unique(d[,repColnm])]
      ### Merge with data:
      d <- rbind(d,dd)
      ### Give a warning:
      warning(paste0('Observation-ID ', levels(x[,IDvar])[i], ' did not have rows for all levels of ', repColnm, '.'))
    }
    ### Order according to repCol:
    d <- d[order(d[,repColnm]),]
    ### Create the wide-format data:
    dwf0 <- d[, !(colnames(d) %in% repVars) & colnames(d)!=repColnm & colnames(d)!=IDvar, drop=FALSE]   # Only fixed variables
    dwf <- cbind(d[,IDvar, drop=FALSE], dwf0)[1,,drop=FALSE]   # Put IDvar in beginning
    dwr <- d[, (colnames(d) %in% repVars), drop=FALSE]   # Only repeated variables
    ### Iterate through the repVar columns:
    djL <- list()
    for(j in 1:ncol(dwr)){
      ### wide data format:
      dj <- do.call(data.frame, as.list(dwr[,j]))
      ### Add appropriate colnames:
      cn0 <- expand.grid(colnames(dwr)[j], lvs)
      cn1 <- do.call(paste0, list(cn0[,2], colNmStr, cn0[,1]))
      colnames(dj) <- cn1
      djL[[j]] <- dj
    }
    ### Combine data:
    dw <- cbind(dwf, do.call(cbind, djL))
    xwL[[i]] <- dw
  }
  ### Merge into final data frame:
  xw <- do.call(rbind, xwL)
  return(xw)
}


#*********************************************************************************
#   MAKE COLOR TRANSPARENT    ####
#*********************************************************************************
mktransp <- function(color, alpha=100){
  stopifnot(alpha >= 0 & alpha <= 255)
  newCol<-col2rgb(color)
  transCol <- apply(newCol, 2, function(x){rgb(red=x[1], green=x[2],
                                               blue=x[3], alpha=alpha,
                                               maxColorValue=255)})
  return(transCol)
}


#*********************************************************************************
#   ADD IMAGE TO PLOT    ####
#*********************************************************************************
addImgToPlot <- function(filename, x=NULL, y = NULL, width = NULL, angle=0, interpolate = TRUE){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  ### Turn image to array object:
  if(grepl(pattern = '.png', filename)){
    obj <- png::readPNG(source = filename)
  }else if(grepl(pattern = '.jpg', filename)){
    obj <- jpeg::readJPEG(source = filename)
  }else{
    stop('You must provide a filepath ending with .png or .jpg')
  }
  ### Get some of the relevant measures:
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  height <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  #' Have to adapt position of center in case of rotation
  #' because rotation is per default around bottom left corner.
  ### Angle of center relative to left bottom before rotation:
  ang0 <- atan(HEIi/WIDi)*180/pi
  ### Distance from corner to center in inches:
  z <- sqrt((HEIi/2)^2 + (WIDi/2)^2)
  ### Angle after rotation
  ang1 <- ang0 + angle
  ### Coordinates of center relative to bottom-left corner after rotation:
  ### In inches:
  xaR_inch <- z * cos(ang1 * pi /180)
  yaR_inch <- z * sin(ang1 * pi /180)
  ### Convert to units:
  xaR <- xaR_inch/PIN[1]*(USR[2]-USR[1])
  yaR <- yaR_inch/PIN[2]*(USR[4]-USR[3])
  ### Coordinates of center relative to bottom-left corner before rotation:
  xbR <- width/2
  ybR <- height/2
  ### Shift of center through rotation:
  xshft <- xbR - xaR
  yshft <- ybR - yaR
  ### Draw image:
  rasterImage(image = obj,
              xleft = (x+xshft)-(width/2), xright = (x+xshft)+(width/2),
              ybottom = (y+yshft)-(height/2), ytop = (y+yshft)+(height/2),
              interpolate = interpolate, angle = angle)
}


#*********************************************************************************
#   PREVALENCE TABLE    ####
#*********************************************************************************
prevTabl <- function(X, FUN, catVar=NULL, atLeastOnce=FALSE){
  ### Create default catVar:
  if(is.null(catVar)){
    catVar <- factor(rep('', nrow(X)))
  }
  stopifnot(is.factor(catVar))   # Must be a factor
  ### Split data according to groups:
  Xsp <- split(X, f = catVar)
  ### Apply function to data:
  tfTbl <- lapply(Xsp, FUN = FUN)
  ### Count successes:
  xs <- lapply(tfTbl, FUN = function(x){apply(x, 2, sum, na.rm=TRUE)})
  xstbl <- do.call(rbind, xs)
  ### Number of entries:
  xn <- lapply(tfTbl, FUN=function(x){apply(x, 2, function(y){sum(!is.na(y), na.rm = TRUE)})})
  xntbl <- do.call(rbind, xn)
  ### Calculate the rates:
  maprt <- Map(f = function(xx,nn){xx/nn}, xx=xs, nn=xn)
  xrat <- do.call(rbind, maprt)   # rbind results
  ### Calculate "onceTrue" column:
  ### Count successes:
  xonce_tf <- lapply(tfTbl, FUN = function(x){apply(x,1, any, na.rm=TRUE)})
  xonce_s <- lapply(xonce_tf, sum, na.rm=TRUE)
  xonce_stbl <- do.call(rbind, xonce_s)
  ### Number of entries:
  xonce_n <- lapply(xonce_tf, function(x){sum(!is.na(x))})
  xonce_ntbl <- do.call(rbind, xonce_n)
  ### Calculate rate:
  xonce_maprt <- Map(f = function(xx,nn){xx/nn}, xx=xonce_s, nn=xonce_n)
  xonce_rat <- do.call(rbind, xonce_maprt)
  colnames(xonce_rat) <- 'atLeastOnce'
  ### Create final table:
  mapfin <- Map(f = function(x,n){paste0(x, '/', n, ' (', round(x/n*100, 1),'%)')}, x=xs, n=xn)
  fintbl <- noquote(do.call(rbind, mapfin))
  colnames(fintbl) <- colnames(X)
  ### Add "onceTrue" column:
  if(atLeastOnce){
    ### Generate final output:
    maponcefin <- Map(f = function(x,n){paste0(x, '/', n, ' (', round(x/n*100, 1),'%)')}, x=xonce_s, n=xonce_n)
    oncefintbl <- noquote(do.call(rbind, maponcefin))
    colnames(oncefintbl) <- colnames(xonce_rat)
    ### Add to final table:
    fintbl <- noquote(cbind(fintbl, oncefintbl))
  }
  ### Collect output:
  if(atLeastOnce){
    L <- list(xs=xstbl, xn=xntbl, xrat=xrat, xonce_s=xonce_stbl,
              xonce_n=xonce_ntbl, xonce_rat=xonce_rat, fintbl=fintbl)
  }else{
    L <- list(xs=xstbl, xn=xntbl, xrat=xrat, fintbl=fintbl)
  }
  ### Print final table:
  print(fintbl)
  ### Silent return of result-list:
  invisible(L)
}


#*********************************************************************************
#   STANDARDIZE DATA TABLE    ####
#*********************************************************************************
standizDat <- function(x, chngName = TRUE, onlyCenter = FALSE){
  ### Change colnames of data:
  if(chngName){
    tag <- ifelse(onlyCenter, '_cntr', '_stnd')
    colnames(x)[sapply(x, is.numeric)] <- paste0(colnames(x)[sapply(x, is.numeric)],
                                                 tag)
  }
  ### Function to be used in lapply:
  lapp_f <- function(y){
    if(is.numeric(y)){
      rval <- scale(y, scale = !onlyCenter)
    } else{
      rval <- y
    }
    return(rval)
  }
  ### Lapply call:
  lap_res <- lapply(x, lapp_f)
  ### combine to dataframe:
  rval <- do.call(data.frame, lap_res)
  return(rval)
}


#*********************************************************************************
#   REMOVE OUTLIERS FROM DATA    ####
#*********************************************************************************
rmOutliers <- function(x, chngName = TRUE){
  ### Change colnames of data:
  if(chngName){
    colnames(x)[sapply(x, is.numeric)] <- paste0(colnames(x)[sapply(x, is.numeric)],
                                                 '_noOutl')
  }
  ### Function to be used in lapply:
  lapp_f <- function(y){
    if(is.numeric(y)){
      ### Find outliers with boxplot func:
      outl <- boxplot.stats(y)$out
      outl_id <- which(y %in% outl)
      y[outl_id] <- NA
      rval <- y
    } else{
      rval <- y
    }
    return(rval)
  }
  ### Lapply call:
  lap_res <- lapply(x, lapp_f)
  ### combine to dataframe:
  rval <- do.call(data.frame, lap_res)
  return(rval)
}


#*********************************************************************************
#   PROGRESS BAR FOR LOOP    ####
#*********************************************************************************
checkprogress <- function(val, endi, starti = 1){
  setTxtProgressBar(txtProgressBar(min=starti-1, max=endi, style = 3), value = val)
}


#*********************************************************************************
#   FLATTEN RECURSIVE LIST    ####
#*********************************************************************************
#' Taken from internet
#' (https://stackoverflow.com/questions/47603578/flatten-recursive-list)
flatten <- function(x) {
  if (!inherits(x, "list")) return(list(x))
  else return(unlist(c(lapply(x, flatten)), recursive = FALSE))
}
