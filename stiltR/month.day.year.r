month.day.year<-function(jul, origin.){
#returns month, day, year
#
#  $Id: month.day.year.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

   if(missing(origin.) || is.null(origin.))
           if(is.null(origin. <- .Options$chron.origin))
                   origin. <- c(month = 1, day = 1, year = 1960)
   if(all(origin. == 0))
           shift <- 0
   else shift <- julian(origin = origin.)
   # relative origin
   # "absolute" origin
   j <- jul + shift
   j <- j - 1721119
   y <- (4 * j - 1) %/% 146097
   j <- 4 * j - 1 - 146097 * y
   d <- j %/% 4
   j <- (4 * d + 3) %/% 1461
   d <- 4 * d + 3 - 1461 * j
   d <- (d + 4) %/% 4
   m <- (5 * d - 3) %/% 153
   d <- 5 * d - 3 - 153 * m
   d <- (d + 5) %/% 5
   y <- 100 * y + j
   y <- y + ifelse(m < 10, 0, 1)
   m <- m + ifelse(m < 10, 3, -9)
   list(month = m, day = d, year = y)
}
