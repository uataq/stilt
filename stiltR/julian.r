julian<-function(m, d, y, origin.){
#returns day since 1/1/1960
#
#  $Id: julian.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

   only.origin <- all(missing(m), missing(d), missing(y))
   if(only.origin)
           m <- d <- y <- NULL
   # return days since origin
   nonnumeric.p <- !c(is.numeric(m), is.numeric(d), is.numeric(y))
   if(any(nonnumeric.p) && !only.origin) {
           badarg <- paste(c("m", "d", "y")[nonnumeric.p], collapse = ", "
                   )
           stop(paste("Arguments:", badarg, "are not numeric"))
   }
   if(missing(origin.))
           if(is.null(origin. <- .Options$chron.origin))
                   origin. <- c(month = 1, day = 1, year = 1960)
   nms <- names(d)
   max.len <- max(length(m), length(d), length(y))
   #
   # prepend new origin value and rep out to common max. length:
   m <- c(origin.[1], rep(m, length = max.len))
   d <- c(origin.[2], rep(d, length = max.len))
   y <- c(origin.[3], rep(y, length = max.len))
   #
   # code from julian date in the S book (p.269)
   #
   y <- y + ifelse(m > 2, 0, -1)
   m <- m + ifelse(m > 2, -3, 9)
   c <- y %/% 100
   ya <- y - 100 * c
   out <- (146097 * c) %/% 4 + (1461 * ya) %/% 4 + (153 * m + 2) %/% 5 +
           d + 1721119
   #
   # now subtract the new origin from all dates
   #
   if(!only.origin) {
           if(all(origin. == 0))
                   out <- out[-1]
           else out <- out[-1] - out[1]
   }
   names(out) <- nms
   out
}
