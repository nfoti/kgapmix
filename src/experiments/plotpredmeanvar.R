plotpredmeanvar <- function(matfn, nstd)
{
  require(ggplot2)
  require(reshape)
  require(R.matlab)
  require(RColorBrewer)

  # Read mat file and parse data
  mdat <- readMat(matfn)

  mnames <- names(mdat)
  xstr <- mdat$xstr[1] # Peculiarity of readMat
  ystr <- mdat$ystr[1]
  N <- length(mdat[['times']])

  times <- mdat[['times']]

  dat <- data.frame(times=times)
  dat.std.min <- data.frame(times=times)
  dat.std.max <- data.frame(times=times)

  mfn <- NULL

  for(n in mnames)
  {
    if(substr(n,1,3) == 'mf.') {
      S <- strsplit(n, "\\.")[[1]]
      if(length(S) != 2)
      {
        cat('Error: Got something weird:',S,'\n',sep=' ')
      }
      dat$newcol <- as.vector(mdat[[eval(n)]])
      names(dat)[length(names(dat))] <- S[2]
      mfn <- c(mfn, S[2])
    } else if(substr(n,1,3) == 'mfv') {
      S <- strsplit(n, "\\.")[[1]]
      if(length(S) != 2)
      {
        cat('Error: Got something weird:',S,'\n',sep=' ')
      }
      dat.std.min$newcol <- as.vector(mdat[[eval(n)]])
      names(dat.std.min)[length(names(dat.std.min))] <- S[2]
      dat.std.max$newcol <- as.vector(mdat[[eval(n)]])
      names(dat.std.max)[length(names(dat.std.max))] <- S[2]
      mfn <- c(mfn, S[2])
    }
  }
  mfn <- unique(mfn)

  for(n in mfn)
  {
    dat.std.min[[eval(n)]] <- dat[[eval(n)]] - nstd*sqrt(dat.std.min[[eval(n)]])
    dat.std.max[[eval(n)]] <- dat[[eval(n)]] + nstd*sqrt(dat.std.max[[eval(n)]])
  }

  # Melt data frame for nice plotting
  dat.melt <- melt(dat, id="times")
  dat.std.min.melt <- melt(dat.std.min, id="times")
  dat.std.max.melt <- melt(dat.std.max, id="times")
  dat.std.melt <- cbind(dat.std.min.melt, dat.std.max.melt$value)
  names(dat.std.melt)[length(names(dat.std.melt))] <- 'vmax'
  names(dat.std.melt)[length(names(dat.std.melt))-1] <- 'vmin'
  names(dat.melt)[names(dat.melt)=='variable'] <- 'sampler'
  names(dat.std.melt)[names(dat.std.melt)=='variable'] <- 'sampler'

  levels(dat.melt$sampler) <- sort(levels(dat.melt$sampler))
  levels(dat.std.melt$sampler) <- sort(levels(dat.std.melt$sampler))
  #if(length(levels(dat.melt$sampler)) == 2) {
  #  levels(dat.melt$sampler) <- c(levels(dat.melt$sampler),'sngp')
  #}
  #if(length(levels(dat.std.melt$sampler)) == 2) {
  #  levels(dat.std.melt$sampler) <- c(levels(dat.std.melt$sampler),'sngp')
  #}
  #cat(length(levels(dat.melt$sampler)),'\n')
  #cat(length(levels(dat.std.melt$sampler)),'\n')

  # Need a separate data frame for the points
  pdat <- data.frame(x=mdat[['X']], y=mdat[['Y']])
  xmin <- min(times)
  xmax <- max(times)
  pdat <- subset(pdat, x >= xmin & x <= xmax)

  mycols <- brewer.pal(3,"Set1")
  names(mycols) <- levels(dat.melt$sampler)

  # Now plot
  p <- ggplot()
  p <- p + geom_point(data=pdat,aes(x=x,y=y),alpha=I(0.2))
  p <- p + geom_line(data=dat.melt,aes(x=times,y=value,colour=sampler),
                     size=I(1))
  p <- p + geom_ribbon(data=dat.std.melt,
                       aes(x=times,ymin=vmin,ymax=vmax,colour=sampler),
                       alpha=I(0.1), size=I(.5), linetype=0)
  p <- p + geom_ribbon(data=dat.std.melt,
                       aes(x=times,ymin=vmin,ymax=vmax,colour=sampler),
                       alpha=I(0.0), size=I(.5), show_guide=FALSE)
  p <- p + theme_bw()
  p <- p + scale_x_continuous(name=xstr,limits=c(xmin,xmax)) + scale_y_continuous(ystr)
  p <- p + scale_colour_manual(name="sampler", values=mycols)

  # This looked good on the screen but not at the scale for the paper
  #p <- p + geom_point(data=pdat,aes(x=x,y=y),alpha=I(0.5))
  #p <- p + geom_line(data=dat.melt,aes(x=times,y=value,colour=sampler),
  #                   size=I(2))
  #p <- p + geom_ribbon(data=dat.std.melt,
  #                     aes(x=times,ymin=vmin,ymax=vmax,colour=sampler),
  #                     alpha=I(0.15), size=I(1.1), linetype=0)
  #p <- p + geom_ribbon(data=dat.std.melt,
  #                     aes(x=times,ymin=vmin,ymax=vmax,colour=sampler),
  #                     alpha=I(0.0), size=I(1.1), show_guide=FALSE)
  #p <- p + theme_bw()
  #p <- p + scale_x_continuous(name=xstr,limits=c(xmin,xmax)) + scale_y_continuous(ystr)
  #p <- p + scale_colour_manual(name="sampler", values=mycols)


  # Didn't work very well
  #p <- p + opts(axis.text.x=theme_text(size=14))
  #p <- p + opts(axis.title.x=theme_text(size=16))
  #p <- p + opts(axis.text.y=theme_text(size=14))
  #p <- p + opts(axis.title.y=theme_text(size=16,vangle=90,vjust=0.0))
  show(p)

  return(p)

}
