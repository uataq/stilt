#
#  $Id: image.plot.fix.r,v 1.5 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

image.plot.fix<-function (..., add = FALSE, nlevel = 64, legend.shrink = 0.9,
    legend.width = 0.05, graphics.reset = FALSE, horizontal = FALSE, 
    offset = 1 * legend.width, bigplot = NULL, smallplot = NULL, 
    legend.only = FALSE, col = topo.colors(nlevel),lg=F) 
{
    old.par <- par(no.readonly = TRUE)
 if(!("package:fields"%in%search())){library(fields);print("loaded library fields")}
    info <- image.plot.info(...)
    if (add) 
        big.plot <- old.par$plt
    if (legend.only) 
        graphics.reset <- TRUE
    temp <- image.plot.plt.fix(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, horizontal = horizontal, 
        offset = offset, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        image(..., add = add, col = col)
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    temp <- list(...)
    iy <- seq(info$zlim[1], info$zlim[2], , nlevel)
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    ix <- 1
    if (!horizontal) {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
#        axis(4, mgp = c(3, 0, 0))
        if(lg){axis(4,labels=paste("1E",axTicks(4)),at=axTicks(4))
          }else{axis(4)}
    }
    else {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
            col = col)
        if(lg){axis(1,labels=paste("1E",axTicks(1)))
          }else{axis(1)}

    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}
