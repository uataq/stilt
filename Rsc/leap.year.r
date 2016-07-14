#
#  $Id: leap.year.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

leap.year<-function(y){
        if(inherits(y, "dates"))
                y <- month.day.year(as.numeric(y), origin = origin(y))$year
        y %% 4 == 0 & (y %% 100 != 0 | y %% 400 == 0)
}
