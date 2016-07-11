day.of.week<-function(month, day, year){
#returns day of week as number (0: Sun, 6: Sat)
#
#  $Id: day.of.week.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

        ix <- year + trunc((month - 14)/12)
        jx <- trunc((13 * (month + 10 - (month + 10) %/% 13 * 12) - 1)/5) +
                day + 77 + (5 * (ix - (ix %/% 100) * 100)) %/% 4 + ix %/% 400 -
                (ix %/% 100) * 2
        jx %% 7
}
