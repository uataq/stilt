read.sgp <- function(fname,co2.name=c('co2.ppm'),range.co2=c(0,2000),site.name='site') {
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------
  if (file.exists(fname)) {
    indat <- as.matrix(read.table(fname,head=TRUE))
    len.all <- dim(indat)[1]
    match.site <- match(site.name,dimnames(indat)[[2]],nomatch=0)
    indat <- indat[,-match.site,drop=F]
    dim.save <- dim(indat)
    names.save <- dimnames(indat)
    indat <- array(as.numeric(indat),dim=dim.save,dimnames=names.save)
  } else {
    len.all <- 0
  }
  sel <- F
  if (len.all > 0) {
    sel <- !is.na(indat[,co2.name])
    if (any(!is.logical(range.co2)))   {
      sel <- sel & (indat[,co2.name] >= range.co2[1])
      sel <- sel & (indat[,co2.name] <= range.co2[2])
    }
  }
  if (sum(sel) == 0) {
    indat <- NULL
  } else {
    indat <- indat[sel,]
  }
  list(data=indat,fname=fname,range.co2=range.co2,len.all=len.all)
}

read.wlef <- function(fname,co2.name=c('co2'),range.co2=c(0,2000)) {
# read data from Argyle file
# optionally subselect data that satisfies criteria of:
#   range.co2: allowable range for mean co2
#
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------
  if (file.exists(fname)) {
    indat <- as.matrix(read.table(fname,head=TRUE))
    len.all <- dim(indat)[1]
  } else {
    len.all <- 0
  }
  sel <- F
  if (len.all > 0) {
    sel <- !is.na(indat[,co2.name])
    if (any(!is.logical(range.co2)))   {
      sel <- sel & (indat[,co2.name] >= range.co2[1])
      sel <- sel & (indat[,co2.name] <= range.co2[2])
    }
  }
  if (sum(sel) == 0) {
    indat <- NULL
  } else {
    indat <- indat[sel,]
  }
  list(data=indat,fname=fname,range.co2=range.co2,len.all=len.all)
}

read.argyle <- function(fname,col.names=c('year','month','day','hour','minute','second',
                             'level','co2','sdco2.tot','sdco2.anal','qcflag'),
                     okflags=c(999,910,900),range.co2=c(0,2000),range.sdco2=c(0,2000),level=F) {
# read data from Argyle file
# optionally subselect data that satisfies criteria of:
#   okflags: acceptable qcflags
#   range.co2: allowable range for mean co2 (default catches -999 values)
#   range.sdco2: allowable range for total uncertainty (default catches -999 values)
#   level: level of tower (1-12m; 2-107m)
#
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------
  if (file.exists(fname)) {
    indat <- scan(fname)
    len.all <- length(indat)/length(col.names)
  } else {
    len.all <- 0
  }
  sel <- F
  if (len.all > 0) {
    dim(indat) <- c(length(col.names),len.all)
    indat <- t(indat)
    dimnames(indat) <- list(NULL,col.names)
    sel <- rep(T,dim(indat)[1])
    if (!is.logical(level)) sel <- sel & indat[,'level'] == level
    if (any(!is.logical(range.co2)))   {
      sel <- sel & (indat[,'co2'] >= range.co2[1])
      sel <- sel & (indat[,'co2'] <= range.co2[2])
    }
    if (any(!is.logical(range.sdco2))) {
      sel <- sel & indat[,'sdco2.tot'] >= range.sdco2[1] 
      sel <- sel & indat[,'sdco2.tot'] <= range.sdco2[2]
    }
    sel.qc <- rep(F,dim(indat)[1])
    do.selqc <- F
    for (tmp.flag in okflags) {
      if (!is.logical(tmp.flag)) {
        do.selqc <- T
        sel.qc <- sel.qc | indat[,'qcflag'] == tmp.flag
      }
    }
    if (do.selqc) sel <- sel & sel.qc
  }
  if (sum(sel) == 0) {
    indat <- NULL
  } else {
    indat <- indat[sel,]
  }
  list(data=indat,fname=fname,okflags=okflags,range.co2=range.co2,range.sdco2=range.sdco2,level=level,len.all=len.all)
}

match.co2 <- function(sim.data,data.list=list(),avg.sec=3600,
                      fname.root='/project/p1229/Mission2004/Data/ARGYLE/FINALDATA/CO2/2004/co2.meancal.',
                      fname.wlef1='/project/p1229/jel/WLEF/wlef',
                      fname.wlef2='hrlyco2.tbl',
                      fname.sgp1='/anvil/carbon1/SGP/sgp',
                      fname.sgp2='.co2.',
                      fname.sgp3='m.dat',
                      start=NULL,
                      eps.latlon=0.02,eps.agl=2,full.list=F,timing=F,...) {
  #function to compute average of measured CO2 that matches the computed
  #needed input: sim.data: output from Trajecvprm
  #                  needed here: columns ident, latstart, lonstart, aglstart
  #optional input:
  # start - matrix with lat,lon,agl that must match receptor locations, and level
  #         that determines which level to read from obs data
  #         dimnames of rows determines which read routine to use
  #         so far only supporting read.argyle and read.wlef
  # ...   - passed to read routine(s)
#
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------
  ptime <- 0
  if (is.null(start)) {
    start <- cbind(lat=45.03,lon=-68.68,agl=10,level=1)
    dimnames(start)[[1]] <- "argyle.1"
    start <- rbind(start,
                   argyle.2=c(45.03,-68.68,107,2))
    start <- rbind(start,
                   wlef.6=c(45.9549,-90.2723,396,1))
    start <- rbind(start,
                   wlef.5=c(45.9549,-90.2723,244,1))
    start <- rbind(start,
                   wlef.4=c(45.9549,-90.2723,122,1))
    start <- rbind(start,
                   wlef.3=c(45.9549,-90.2723,076,1))
    start <- rbind(start,
                   wlef.2=c(45.9549,-90.2723,030,1))
    start <- rbind(start,
                   wlef.1=c(45.9549,-90.2723,011,1))
    start <- rbind(start,
                   sgp.1=c(36.61,-97.49,60,1))
  }
  co2 <- rep(NA,dim(sim.data)[1])
  nco2 <- rep(0,dim(sim.data)[1])
  dnames1 <- rep('',dim(sim.data)[1])
  dnames2 <- rep('',dim(sim.data)[1])
  if (timing) ptime <- proc.time()
  for (iday in 1:dim(sim.data)[1]) {
    if (any(is.na(sim.data[iday,c('latstart','lonstart','aglstart')]))) {
      sel.start <- F
    } else {
      sel.start <- (abs(sim.data[iday,'latstart']-start[,'lat']) <= eps.latlon &
                    abs(sim.data[iday,'lonstart']-start[,'lon']) <= eps.latlon &
                    abs(sim.data[iday,'aglstart']-start[,'agl']) <= eps.agl)
    }
    if (sum(sel.start) == 1) {
      ident <- sim.data[iday,'ident']
      jul.1 <- ident-0.5*avg.sec/86400
      jul.2 <- ident+0.5*avg.sec/86400
      ymd.1 <- month.day.year(floor(jul.1))
      ymd.2 <- month.day.year(floor(jul.2))
      yyyymmdd.1 <- paste(sprintf('%4d',as.integer(ymd.1$year)),
                          sprintf('%2.2d',as.integer(ymd.1$month)),sprintf('%2.2d',as.integer(ymd.1$day)),sep='')
      yyyymmdd.2 <- paste(sprintf('%4d',as.integer(ymd.2$year)),
                          sprintf('%2.2d',as.integer(ymd.2$month)),sprintf('%2.2d',as.integer(ymd.2$day)),sep='')
      yyyymmdd <- unique(c(yyyymmdd.1,yyyymmdd.2))
      measured <- 0
      nmeasured <- 0
      i.tmp <- 0
      dname.previous <- ''
      tmp.dat <- NULL
      tmp.jul <- NULL
      if (!is.na(pmatch('argyle',dimnames(start)[[1]][sel.start]))) {
        #for argyle data:
        for (tmp.name in yyyymmdd) {
          i.tmp <- i.tmp + 1
          dname <- paste(dimnames(start)[[1]][sel.start],tmp.name,sep='.')
                                        # read data if needed
          if (i.tmp == 1) dnames1[iday] <- dname
          if (i.tmp == 2) dnames2[iday] <- dname
          if (is.null(data.list[[dname]])) {
            if (!is.na(pmatch('argyle',dname))) {
              data.list[[dname]] <- read.argyle(fname=paste(fname.root,tmp.name,sep=''),
                                                level=start[sel.start,'level'],...)
              if (timing) cat('Timing: read.argyle=',proc.time()-ptime,'\n')
            }
          }
          tmp.dat <- data.list[[dname]]$data
          if (!is.null(tmp.dat)) {
            tmp.jul <- julian(y=tmp.dat[,'year'],m=tmp.dat[,'month'],d=tmp.dat[,'day']) +
              (tmp.dat[,'hour'] + (tmp.dat[,'minute']+tmp.dat[,'second']/60.)/60.)/24.
            sel.avg <- tmp.jul >= jul.1 & tmp.jul <= jul.2 #if all F, zeroes are added here:
            measured <- measured + sum(tmp.dat[sel.avg,'co2'])
            nmeasured <- nmeasured + sum(sel.avg)
          }
        } #endfor tmp.name
      } #endif argyle data
      if (!is.na(pmatch('wlef',dimnames(start)[[1]][sel.start]))) {
        #for WLEF data:
        yyyy <- unique(substring(yyyymmdd,1,4))
        for (tmp.name in yyyy) {
          i.tmp <- i.tmp + 1
          dname <- paste(dimnames(start)[[1]][sel.start],tmp.name,sep='.')
                                        # read data if needed
          if (i.tmp == 1) dnames1[iday] <- dname
          if (i.tmp == 2) dnames2[iday] <- dname
          if (dname != dname.previous) {
            if (is.null(data.list[[dname]])) {
              if (!is.na(pmatch('wlef',dname))) {
                data.list[[dname]] <- read.wlef(fname=paste(fname.wlef1,
                                                sprintf("%3.3i",start[sel.start,'agl']),'.',
                                                tmp.name,fname.wlef2,sep=''),...)
                if (timing) cat('Timing: read.wlef=',proc.time()-ptime,'\n')
              }
            }
            tmp.dat <- data.list[[dname]]$data
            tmp.jul <- NULL
            if (!is.null(tmp.dat)) tmp.jul <-
              julian(y=tmp.dat[,'year'],m=tmp.dat[,'month'],d=tmp.dat[,'day']) +
                (tmp.dat[,'hr'])/24.
            dname.previous <- dname
          }
          if (!is.null(tmp.jul)) {
            sel.avg <- tmp.jul >= jul.1 & tmp.jul <= jul.2 #if all F, zeroes are added here:
            measured <- measured + sum(tmp.dat[sel.avg,'co2'])
            nmeasured <- nmeasured + sum(sel.avg)
          }
        } #endfor tmp.name
      } #endif WLEF data
      if (!is.na(pmatch('sgp',dimnames(start)[[1]][sel.start]))) {
        #for ARM SGP data:
        yyyy <- unique(substring(yyyymmdd,1,4))
        for (tmp.name in yyyy) {
          i.tmp <- i.tmp + 1
          dname <- paste(dimnames(start)[[1]][sel.start],tmp.name,sep='.')
                                        # read data if needed
          if (i.tmp == 1) dnames1[iday] <- dname
          if (i.tmp == 2) dnames2[iday] <- dname
          if (dname != dname.previous) {
            if (is.null(data.list[[dname]])) {
              if (!is.na(pmatch('sgp',dname))) {
                data.list[[dname]] <- read.sgp(fname=paste(fname.sgp1,tmp.name,fname.sgp2,
                                                sprintf("%2.2i",start[sel.start,'agl']),
                                                fname.sgp3,sep=''),co2.name='co2.ppm',...)
                if (timing) cat('Timing: read.sgp=',proc.time()-ptime,'\n')
              }
            }
            tmp.dat <- data.list[[dname]]$data
            tmp.jul <- NULL
            if (!is.null(tmp.dat)) tmp.jul <-
              julian(y=tmp.dat[,'yr'],m=tmp.dat[,'mo'],d=tmp.dat[,'day']) +
                (tmp.dat[,'hr'])/24.
            dname.previous <- dname
          }
          if (!is.null(tmp.jul)) {
            sel.avg <- tmp.jul >= jul.1 & tmp.jul <= jul.2 #if all F, zeroes are added here:
            measured <- measured + sum(tmp.dat[sel.avg,'co2.ppm'])
            nmeasured <- nmeasured + sum(sel.avg)
          }
        } #endfor tmp.name
      } #endif sgp data
      if (nmeasured > 0) {
        measured <- measured/nmeasured
      } else {
        measured <- NA
      }
      co2[iday] <- measured
      nco2[iday] <- nmeasured
    } #endif sel.start
    if (timing) {
      ntime <- proc.time()
      cat('Timing: iday',iday,'tmp.jul:',length(tmp.jul),'=',ntime-ptime,'\n')
      ptime <- ntime
    }
  }
#return average CO2 measurements, and updated data.list
  if (full.list) {
    out.list <- list(co2=co2,nco2=nco2,data.list=data.list,dnames1=dnames1,dnames2=dnames2)
  } else {
    out.list <- list(co2=co2,nco2=nco2,dnames1=dnames1,dnames2=dnames2)
  }
  out.list
}

plot.co2 <- function(sim.data,match.data=NULL,gee.names=NULL,resp.names=NULL,
                     bgname='co2ini',ffname='co2ffm',
                     add=F,axis.julian=T,axis.time=T,ylim=NULL,xlim=NULL,ylab='CO2 [ppm]',
                     type='b',pch=1:2,col=1:2,lty=1:2,do.plot=T,format='%m/%d/%H',xlab='Time',
                     at=NULL,at.julian=at,at.time=NULL,by.time="12 hour",cex=1,lwd=1,...) {
  # plot up simulated and observed CO2, and returned matched up vectors of each
  # Required input:
  # sim.data  - output from Trajecvprm (2d array)
  #             if more than one object is to be plotted, provide a character vector
  #             of R object names (in search path)
  # Optional inputs:
  #    only if add=F: match.data: observed data matched to sim.data
  #                   if not supplied, will be generated by a call to match.co2
  #    gee.names, resp.names, bgname, ffname: column names for gee, respiration,
  #                   background, and fossil fuel contributions to simulated CO2
  #    axis.time - (T/F) use R date object for plotting x axis
  #    format - format string (see strptime) for time axis plotting
  #    ... - passed on to plot routine
#
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------
  if (is.character(sim.data)) {
    tmp.sim <- sim.data
    sim.data <- NULL
    for (x in tmp.sim) sim.data <- rbind(sim.data,get(x))
  }

  if(is.null(gee.names) || is.null(resp.names)) {
    gee.names <- character(0);resp.names <- character(0)
    for (x in dimnames(sim.data)[[2]]) {
      if (!is.na(pmatch('gee',x))) gee.names <- c(gee.names,x)
      if (!is.na(pmatch('resp',x))) resp.names <- c(resp.names,x)
    }
  }
  if (length(gee.names) != length(resp.names)) warning(paste('unequal lengths of gee.names=',
              paste(gee.names,collapse=';'),'and resp.names=',paste(resp.names,collapse=';')))
  sim.co2 <- rep(0,dim(sim.data)[1])
  if (length(bgname) > 0) sim.co2 <- sim.co2 + sim.data[,bgname]
  if (length(ffname) > 0) sim.co2 <- sim.co2 + sim.data[,ffname]
  if (length(gee.names) > 0) sim.co2 <- sim.co2 + apply(sim.data[,gee.names,drop=F],MAR=1,FUN=sum)
  if (length(resp.names) > 0) sim.co2 <- sim.co2 + apply(sim.data[,resp.names,drop=F],MAR=1,FUN=sum)
  if (is.null(xlim)) xlim <- range(sim.data[,'ident'])
  xmask <- sim.data[,'ident'] >= xlim[1] & sim.data[,'ident'] <= xlim[2]
  xvals <- sim.data[xmask,'ident']
  xlim.time <- NULL
  xvals.time <- NULL
  if (axis.time) {
    xvals.time <- julian.2.ISO(xvals)
    xlim.time <- julian.2.ISO(xlim)
  }
  obs.co2 <- NULL
  icurve <- 1
  if (!add) {
    if (length(pch) < 2) pch <- rep(pch,2)
    if (length(col) < 2) col <- rep(col,2)
    if (length(lty) < 2) lty <- rep(lty,2)
    if (is.null(match.data)) match.data <- match.co2(sim.data)
    if (is.null(ylim)) ylim <- range(c(sim.co2[xmask],match.data$co2[xmask]),na.rm=TRUE)
    if (do.plot) {
      plot(xlim,ylim,type='n',ylab=ylab,xlab=xlab,xaxt='n',...)
      axis.side <- 1
      if (axis.julian) {
        axis(axis.side,at=at,labels=TRUE)
        axis.side <- 3
      }
      if (axis.time) {
        if (is.null(at.time)) {
          if (is.null(at)) {
            at.time=seq(julian.2.ISO(floor(xlim[1])),xlim.time[2],by=by.time)
            at.time.julian <- ISO.2.julian(at.time)
            at.time <- at.time[at.time.julian >= xlim[1] & at.time.julian <= xlim[2]]
            at.time.julian <- at.time.julian[at.time.julian >= xlim[1] & at.time.julian <= xlim[2]]
          } else {
            at.time.julian <- at
            at.time <- julian.2.ISO(at.time.julian)
          }
        }
        axis(axis.side,at=at.time.julian,labels=format(at.time,format=format,tz="GMT"))        
      }
    }
    obs.co2 <- match.data$co2
    if (do.plot) {
      if (type == 'b' || type == 'p') points(xvals,obs.co2[xmask],pch=pch[icurve],col=col[icurve],cex=cex)
      if (type == 'b' || type == 'l') lines(xvals,obs.co2[xmask],lty=lty[icurve],col=col[icurve],lwd=lwd)
    }
    icurve <- icurve+1
  }
  if (do.plot) {
    if (type == 'b' || type == 'p') points(xvals,sim.co2[xmask],pch=pch[icurve],col=col[icurve],cex=cex)
    if (type == 'b' || type == 'l') lines(xvals,sim.co2[xmask],lty=lty[icurve],col=col[icurve],lwd=lwd)
  }
  invisible(list(xvals=xvals,sim.co2=sim.co2[xmask],obs.co2=obs.co2[xmask],xlim=xlim,
                 xvals.time=xvals.time,xlim.time=xlim.time,ylim=ylim,xmask=xmask))
}

plot.scatter.co2 <- function(in.list=list(),sim.co2=in.list$sim.co2,
                             obs.co2=in.list$obs.co2,
                             gmtrange=NULL,xvals=in.list$xvals,
                             namey='sim',namex='obs',
                             xlab='Observed',ylab='Simulated',
                             add=FALSE,xlim=NULL,ylim=xlim,
                             xleg=xlim[1],yleg=ylim[2],cleg=0.6,
                             tleg=c('N','bias','rms','corr',
                              paste(namey,namex,sep='.'),paste(namex,namey,sep='.')),digits=3,
                             psfile=NULL,cex=0.5,col=1,pch=1,lwd.one=2,col.one=1,reg.lines=TRUE,
                             daily.avg=FALSE,daily.offset=0,...) {
  # plot up scatterplot of simulated and observed CO2
  # Required input:
  # in.list  - output from plot.co2, OR
  # sim.co2  - simulated CO2 (vector)
  # obs.co2  - observed CO2 (vector) 
  #
  # Optional inputs:
  #    gmtrange - only select times within those hours (0-23.99) GMT
  #    xvals    - julian days (needed if gmtrange specified)
  #    daily.avg - flag for computing daily averages before doing scatterplot
  #    daily.offset - by how many hours to offset the Julian day before applying floor()
  #                   For EST (GMT-5 hours), specify offset=-5
#
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------

  xmask <- rep(TRUE,length(sim.co2))
  if (!is.null(gmtrange)) {
    if (is.null(xvals)) stop ('Need to provide xvals if specifying gmtrange\n')
    gmthrs <- xvals %% 1
    if (gmtrange[1] < gmtrange[2]) {
      xmask <- gmthrs >= gmtrange[1]/24. & gmthrs <= gmtrange[2]/24.
    } else {
      xmask <- gmthrs >= gmtrange[1]/24. | gmthrs <= gmtrange[2]/24.
    }
  }
  if (sum(xmask) < 2) {
    cat('Empty selection after applying gmtrange=',gmtrange,' - nothing plotted\n')
    return(list())
  }
  if (daily.avg) {
    if (is.null(gmtrange)) stop('Can only do daily.avg if supplying gmtrange and xvals')
    xdays <- floor(xvals+daily.offset/24)
    gmthrs <- xvals %% 1
    uniqdays <- unique(xdays)
    daily.sim <- rep(NA,length(uniqdays))
    daily.obs <- daily.sim
    nobs <- rep(0,length(uniqdays))
    for (iday in 1:length(uniqdays)) {
      xday <- uniqdays[iday]
      daily.sim[iday] <- mean(sim.co2[xday == xdays & xmask],na.rm=TRUE)
      daily.obs[iday] <- mean(obs.co2[xday == xdays & xmask],na.rm=TRUE)
      nobs[iday] <- sum(!is.na(sim.co2[xday == xdays & xmask]) &
                        !is.na(obs.co2[xday == xdays & xmask]))
    }
    sim.co2 <- daily.sim
    obs.co2 <- daily.obs
  } else {
    nobs <- NULL
    sim.co2 <- sim.co2[xmask]
    obs.co2 <- obs.co2[xmask]
  }
  out.list <- list(sim.co2=sim.co2,obs.co2=obs.co2,psfile=psfile)
  if(!is.null(nobs)) out.list$nobs <- nobs  
  if (!is.null(psfile)) postscript(psfile,paper='special',height=8,width=8,horizontal=FALSE)
  if (is.null(xlim)) xlim <- range(c(sim.co2,obs.co2),na.rm=TRUE)
  if (is.null(ylim)) ylim <- xlim
  out.list$xlim <- xlim
  out.list$ylim <- ylim
  if (!add) {
    plot(obs.co2,sim.co2,xlim=xlim,ylim=xlim,xlab=xlab,ylab=ylab,type='n',...)
    lines(c(max(xlim[1],ylim[1]),min(xlim[2],ylim[2])),c(max(xlim[1],ylim[1]),min(xlim[2],ylim[2])),
          col=col.one,lwd=lwd.one)
  }
  points(obs.co2,sim.co2,cex=cex,col=col,pch=pch)
  sim.on.obs <- lm(sim.co2 ~ obs.co2,list(sim.co2=sim.co2,obs.co2=obs.co2))
  out.list$sim.on.obs=sim.on.obs
  obs.on.sim <- lm(obs.co2 ~ sim.co2,list(sim.co2=sim.co2,obs.co2=obs.co2))
  out.list$obs.on.sim=obs.on.sim
  if (!is.na(xleg) && !is.na(yleg)) {
    l.values <- character(0)
    for (x in tleg) {
      xx <- x
      if (x == paste(namey,namex,sep='.')) xx <- 'sim.obs'
      if (x == paste(namex,namey,sep='.')) xx <- 'obs.sim'
      l.values <- c(l.values,
                    format(switch(xx,
                                  'N'=sum(!is.na(obs.co2) & !is.na(sim.co2)),
                                  'bias'=mean(sim.co2 - obs.co2,na.rm=TRUE),
                                  'rms'=sqrt(mean((sim.co2 - obs.co2)^2,na.rm=TRUE)),
                                  'corr'=cor(sim.co2,obs.co2,use="pairwise.complete.obs"),
                                  'sim.obs'=coef(sim.on.obs)[2],
                                  'obs.sim'=1/coef(obs.on.sim)[2],
                                  NA),
                           digits=digits))
    }
    txtleg <- paste(tleg,'=',l.values)
    legend(xleg,yleg,legend=txtleg,cex=cleg)
    for (x in c('xleg','yleg','txtleg')) out.list[[x]] <- eval(parse(text=x))
  }
  if (reg.lines && all(is.finite(c(coef(sim.on.obs),
                                   -coef(obs.on.sim)[1]/coef(obs.on.sim)[2],
                                   1/coef(obs.on.sim)[2])))) {
    abline(coef(sim.on.obs),lty=1,col=col)
    abline(c(-coef(obs.on.sim)[1]/coef(obs.on.sim)[2],1/coef(obs.on.sim)[2]),lty=2,col=col)
  }
  if(!is.null(psfile)) {
    dev.off()
    cat('Created psfile=',psfile,'\n')
  }
  invisible(out.list)
}

plot.contributions <- function(sim.data,match.data=NULL,gee.names=NULL,resp.names=NULL,
                               bgname='co2ini',ffname='co2ffm',
                               xleg=c(0,0),yleg=c(1,1.0),ncol.leg=c(1,2),cex.leg=c(1,1),
                               xlab='Time',title.top=NULL,title.ts=NULL,title.cont=NULL,cex.title=2,
                               ylim.ts=NULL,ylim.cont=c(-40,40),
                               xlim=NULL,
                               par.list=list(mfrow=c(2,1),xpd=TRUE),
                               do.plot=TRUE,format='%m/%d/%H',
                               at=NULL,at.julian=at,at.time=NULL,
                               by.time="12 hour",by.tick=NULL,tcl.00=-0.7,
                               do.sim=TRUE,do.bg=TRUE,do.ff=TRUE,vtype='net',vtype.leg='',veg=NULL,
                               ylab=expression(paste('CO'[2],' [ppm]')),
                               par.ts=list(),col.ts=1:2,lty.ts=1:2,leg.ts=c('Model','Obs'),
                               col.cont=NULL,pch.cont=NULL,lty.cont=NULL,
                               par.cont=list(),
                               add.ab=c(0.9,0.9),
                               shade.boxes=NULL,shade.density=NULL,shade.col=NA,shade.border=NULL,
                               shade.lty=NULL,shade.lwd=NULL,shade.angle=NULL,
                               arrows.ts=NULL,arrows.cont=NULL,arrows.lwd=2,
                               ...) {
  # plot up contributions to simulated CO2
  # Required input:
  # sim.data  - output from Trajecvprm (2d array)
  # Optional inputs:
  #    gee.names, resp.names, bgname, ffname: column names for gee, respiration,
  #                   background, and fossil fuel contributions to simulated CO2
  #    axis.time - (T/F) use R date object for plotting x axis
  #    format - format string (see strptime) for time axis plotting
  #    ... - passed on to plot routine
#
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------
  if (is.character(sim.data)) stop('Need to provide a matrix for sim.data')
  if (do.plot && !is.null(par.list)) old.par <- par(par.list)

  if(is.null(gee.names) || is.null(resp.names)) {
    gee.names <- character(0);resp.names <- character(0)
    for (x in dimnames(sim.data)[[2]]) {
      if (!is.na(pmatch('gee',x))) gee.names <- c(gee.names,x)
      if (!is.na(pmatch('resp',x))) resp.names <- c(resp.names,x)
    }
  }
  veg.names <- substring(gee.names,4)
  if (is.null(veg)) veg <- veg.names
  veg.names <- veg.names[match(veg,veg.names,nomatch=0)]
  if (length(veg) != length(veg.names)) cat('Some of specified veg=',
              paste(veg,collapse=' '),'omitted, keeping only:',
              paste(veg.names,collapse=' '),'\n')
  veg <- veg.names
  if (length(par.ts) > 0) par.new <- par(par.ts)
  plt.list <- plot.co2(sim.data,match.data,gee.names,resp.names,bgname,ffname,add=FALSE,
                       type='l',lty=lty.ts,col=col.ts,axis.julian=F,axis.time=F,
                       xlim=xlim,ylab=expression(paste("CO"[2]," [ppm]")),xlab='',ylim=ylim.ts,...)
  par.usr <- par('usr')
  if (!is.na(xleg[1]) && !is.na(yleg[1])) {
    legend(par.usr[1] + diff(par.usr[1:2])*xleg[1],par.usr[3] + diff(par.usr[3:4])*yleg[1],
           legend=leg.ts,ncol=ncol.leg[1],col=col.ts,lty=lty.ts,cex=cex.leg[1])
  }
  if (!is.null(arrows.ts)) for (i in 1:length(arrows.ts$x1))
    arrows(par.usr[1]+arrows.ts$x1[i]*diff(par.usr[1:2]),par.usr[3]+arrows.ts$y1[i]*diff(par.usr[3:4]),
	   par.usr[1]+arrows.ts$x2[i]*diff(par.usr[1:2]),par.usr[3]+arrows.ts$y2[i]*diff(par.usr[3:4]),
           lwd=arrows.lwd)

  if(!is.logical(add.ab)) text(par.usr[1]+add.ab[1]*diff(par.usr[1:2]),
                               par.usr[3]+add.ab[2]*diff(par.usr[3:4]),
                               '(a)',cex=cex.title)
  if (is.null(title.ts)) title.ts <- title.top
  if(!is.null(title.ts)) mtext(text=title.ts,side=3,line=2,cex=cex.title)

  if(is.null(xlim)) xlim <- plt.list$xlim
  xlim.time <- julian.2.ISO(xlim)
  at.time=seq(julian.2.ISO(floor(xlim[1])),xlim.time[2],by=by.time)
  at.time.julian <- ISO.2.julian(at.time)
  at.time <- at.time[at.time.julian >= xlim[1] & at.time.julian <= xlim[2]]
  at.time.julian <- at.time.julian[at.time.julian >= xlim[1] & at.time.julian <= xlim[2]]
  
  by.00="24 hour"
  at.time.00=seq(julian.2.ISO(floor(xlim[1])),xlim.time[2],by=by.00)
  at.time.julian.00 <- ISO.2.julian(at.time.00)
  at.time.00 <- at.time.00[at.time.julian.00 >= xlim[1] & at.time.julian.00 <= xlim[2]]
  at.time.julian.00 <- at.time.julian.00[at.time.julian.00 >= xlim[1] & at.time.julian.00 <= xlim[2]]

  axis(side=3,at=at.time.julian,labels=format(at.time,format=format,tz="GMT")) 
  axis(side=3,at=at.time.julian.00,labels=FALSE,tcl=tcl.00)
  if (!is.null(by.tick)) {
    at.time.tick=seq(julian.2.ISO(floor(xlim[1])),xlim.time[2],by=by.tick)
    at.time.julian.tick <- ISO.2.julian(at.time.tick)
    at.time.tick <- at.time.tick[at.time.julian.tick >= xlim[1] & at.time.julian.tick <= xlim[2]]
    at.time.julian.tick <- at.time.julian.tick[at.time.julian.tick >= xlim[1] & at.time.julian.tick <= xlim[2]]
    axis(side=3,at=at.time.julian.tick,labels=FALSE)
  }

  if (!is.null(shade.boxes)) {
    if (is.null(shade.lty)) shade.lty <- par('lty')
    if (is.null(shade.lwd)) shade.lwd <- par('lwd')
    if (is.null(shade.angle)) shade.angle <- 45
    for (i in 1:length(shade.boxes$xleft)) {
      rect(shade.boxes$xleft[i],par.usr[3],shade.boxes$xright[i],par.usr[4],
           density=shade.density,angle=shade.angle,col=shade.col,border=shade.border,
           lty=shade.lty,lwd=shade.lwd)
    }
  }


  if (length(par.ts) > 0) par(par.new)

  if (length(par.cont) > 0) par.new <- par(par.cont)

  yvals <- NULL
  vtypeveg <- paste(vtype.leg,veg,sep='')
  col.names <- c(c('sim-obs','bg-obs','ff')[c(do.sim,do.bg,do.ff)],vtypeveg)
  ncols <- length(col.names)
  if (ncols < 1) stop ('Nothing to be done: do.sim,do.bg,do.ff=F, veg=empty\n')
  yvals <- array(NA,dim=c(sum(plt.list$xmask),ncols),dimnames=list(NULL,col.names))
  if (do.sim) yvals[,'sim-obs'] <- plt.list$sim.co2-plt.list$obs.co2
  if (do.bg) yvals[,'bg-obs'] <- sim.data[plt.list$xmask,bgname]-plt.list$obs.co2
  if (do.ff) yvals[,'ff'] <- sim.data[plt.list$xmask,ffname]
  vtype.match <- FALSE
  if (vtype == 'net') {
    yvals[,vtypeveg] <- sim.data[plt.list$xmask,paste('gee',veg,sep='')] +
      sim.data[plt.list$xmask,paste('resp',veg,sep='')]
    vtype.match=TRUE
  } 
  if (vtype == 'gee') {
    yvals[,vtypeveg] <- sim.data[plt.list$xmask,paste('gee',veg,sep='')]
    vtype.match=TRUE
  } 
  if (vtype == 'resp') {
    yvals[,vtypeveg] <- sim.data[plt.list$xmask,paste('resp',veg,sep='')]
    vtype.match=TRUE
  } 
  if (!vtype.match) stop(paste('invalid vtype=',vtype))
  if(is.null(col.cont)) col.cont <- 1+((1:ncols)-1) %% 6
  if(is.null(pch.cont)) pch.cont <- 2 + ((1:ncols)-1) %/% 6
  if(is.null(lty.cont)) lty.cont <- rep(1,ncols)
  
  if (is.null(ylim.cont)) ylim.cont <- range(yvals,na.rm=TRUE)
  plot(xlim,ylim.cont,type='n',ylab=ylab,xlab=xlab,xaxt='n',...)
  axis(side=1,at=at.time.julian,labels=format(at.time,format=format,tz="GMT")) 
  axis(side=1,at=at.time.julian.00,labels=FALSE,tcl=tcl.00) 
  if (!is.null(by.tick)) axis(side=1,at=at.time.julian.tick,labels=FALSE)
  for (icurve in 1:ncols) {
    points(plt.list$xvals,yvals[,icurve],pch=pch.cont[icurve],col=col.cont[icurve],lty=1)
    lines(plt.list$xvals,yvals[,icurve],lty=lty.cont[icurve],col=col.cont[icurve])
  }
  par.usr <- par('usr')
  if (!is.na(xleg[2])) {
    legend(par.usr[1] + diff(par.usr[1:2])*xleg[2],par.usr[3] + diff(par.usr[3:4])*yleg[2],
           legend=col.names,ncol=ncol.leg[2],pch=pch.cont,col=col.cont,cex=cex.leg[2])
  }
  if(!is.null(title.cont)) mtext(text=title.cont,side=3,line=2,cex=cex.title)
  if(!is.logical(add.ab)) text(par.usr[1]+add.ab[1]*diff(par.usr[1:2]),
                               par.usr[3]+add.ab[2]*diff(par.usr[3:4]),
                               '(b)',cex=cex.title)

  if (!is.null(shade.boxes)) {
    if (is.null(shade.lty)) shade.lty <- par('lty')
    if (is.null(shade.lwd)) shade.lty <- par('lwd')
    if (is.null(shade.angle)) shade.angle <- 45
    for (i in 1:length(shade.boxes$xleft)) {
      rect(shade.boxes$xleft[i],par.usr[3],shade.boxes$xright[i],par.usr[4],
           density=shade.density,angle=shade.angle,col=shade.col,border=shade.border,
           lty=shade.lty,lwd=shade.lwd)
    }
  }
  if (!is.null(arrows.cont)) for (i in 1:length(arrows.cont$x1))
    arrows(par.usr[1]+arrows.cont$x1[i]*cont(par.usr[1:2]),par.usr[3]+arrows.cont$y1[i]*diff(par.usr[3:4]),
	   par.usr[1]+arrows.cont$x2[i]*cont(par.usr[1:2]),par.usr[3]+arrows.cont$y2[i]*diff(par.usr[3:4]),
           lwd=arrows.lwd)

  if (length(par.cont) > 0) par(par.new)
  if (!is.null(par.list)) par(old.par)
}

mean.diurnal.cycle <- function(in.list=list(),jdays=in.list$xvals,tser1=in.list$sim.co2,
                               tser2=in.list$obs.co2,hrstart=0,hrbin=1,jday1=NULL,jday2=NULL) {
  # Function to compute mean diurnal cycle(s) of CO2
  # Input:
  #  Required:
  #    in.list <- list returned by plot.co2
  #  OR
  #    jdays <- Julian days (float)
  #    tser1 <- CO2 values corresponding to jdays
  #    tser2 <- if non-NULL, only use jdays for which both tser1 and tser2 are non-missing
  #  Optional:
  #  hrstart <- hour (GMT) at which to start the diurnal cycle
  #  hrbin <- bin width (hours) of output diurnal cycle
  #  jday1, jday2 <- limits range of jdays to use in computation
#
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------

  if (is.null(jdays) || is.null(tser1)) stop ('Missing required inputs jdays and/or tser1')
  if (length(jdays) != length(tser1)) stop('jdays and tser1 must be of same length')
  if (!is.null(tser2) && length(tser1) != length(tser2)) stop('tser1 and tser2 must be of same length')
  out.list <- list()
  if (all(is.na(jdays)) || all(is.na(tser1))) {
    cat('mean.diurnal.cycle: All missing data, returning empty list\n')
    return(out.list)
  }
  if (is.null(jday1)) jday1 <- min(jdays,na.rm=TRUE)
  if (is.null(jday2)) jday2 <- max(jdays,na.rm=TRUE)
  mask.jdays <- !is.na(jdays)
  mask.jdays <- mask.jdays & jdays >= jday1 & jdays <= jday2
  if(sum(mask.jdays) <= 0) {
    cat(paste('mean.diurnal.cycle: No non-missing dates between specified jday1= ',
              jday1,' and jday2= ',jday2,'\n'))
    return(out.list)
  }
  mask.jdays <- mask.jdays & !is.na(tser1)
  if (!is.null(tser2)) mask.jdays <- mask.jdays & !is.na(tser2)
  if(sum(mask.jdays) <= 0) {
    cat(paste('mean.diurnal.cycle: No non-missing data between specified jday1= ',
              jday1,' and jday2= ',jday2,'\n'))
    return(out.list)
  }
  bin.hr.mids <- seq(hrstart,hrstart+24-hrbin/2,by=hrbin) %% 24
  nbins <- length(bin.hr.mids)
  bin.hr.lims <- c(hrstart-hrbin/2.+(0:nbins)*hrbin) %% 24
  jday.hrs <- (jdays %% 1) * 24
  bin.n <- rep(0,nbins)
  bin.sum1 <- rep(0,nbins)
  bin.sum2 <- rep(0,nbins)
  for (i in 1:nbins) {
    mask.bin <- mask.jdays
    if (bin.hr.lims[i+1] > bin.hr.lims[i]) {
      mask.bin <- (mask.bin &
                   (jday.hrs >= bin.hr.lims[i] & jday.hrs < bin.hr.lims[i+1]))
    } else {
      mask.bin <- (mask.bin &
                   (jday.hrs >= bin.hr.lims[i] | jday.hrs < bin.hr.lims[i+1]))
    }
    bin.n[i] <- bin.n[i] + sum(mask.bin)
    bin.sum1[i] <- bin.sum1[i] + sum(tser1[mask.bin])
    if (!is.null(tser2)) bin.sum2[i] <- bin.sum2[i] + sum(tser2[mask.bin])
  }
  bin.sum1[bin.n <= 0] <- NA
  bin.list <- list(hr=bin.hr.mids,n=bin.n,mean1=bin.sum1/bin.n)
  if (!is.null(tser2)) {
    bin.sum2[bin.n <= 0] <- NA
    bin.list$mean2 <- bin.sum2/bin.n
  }
  out.list$hrstart <- hrstart
  out.list$hrbin <- hrbin
  out.list$jday1 <- jday1
  out.list$jday2 <- jday2
  out.list$bins <- bin.list
  out.list
}

julian.2.ISO <- function(x,...) {
#
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------
  xdmy <- month.day.year(floor(x),...)
  rmdr <- x-floor(x)
  hr <- floor(24*rmdr)
  min <- floor((24*rmdr-hr)*60)
  sec <- round(((24*rmdr-hr)*60-min)*60)
  x.time <- ISOdate(xdmy$year,xdmy$month,xdmy$day,hour=hr,min=min,sec=sec)
  x.time
}

ISO.2.julian <- function(x,...) {
#
#  $Id: measured.co2.r,v 1.6 2008-03-27 14:20:10 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------
  x.char <- format(x,format='%Y_%m_%d_%H_%M_%S',tz="GMT")
  x.julian <- julian(y=as.numeric(substr(x.char,1,4)),
                     m=as.numeric(substr(x.char,6,7)),
                     d=as.numeric(substr(x.char,9,10)),...) +
                       (as.numeric(substr(x.char,12,13)) +
                        (as.numeric(substr(x.char,15,16)) +
                         as.numeric(substr(x.char,18,19))/60.)/60.)/24.

  x.julian
}
