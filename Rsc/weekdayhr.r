weekdayhr<-function(yr,mon,day,hr,runtt,diffGMT=NA){
#Determines 1) day of week   2) hour of day  that would be used in 'Pcogrid' to determine emission factors used to scale CO emissions
#'yr','mon','day','hr' are individual numbers used to specify the starting time
#		let 'hr' run from 0~23
#'runtt' is the runtime from starting time in MINUTES; 'runtt' can be both pos or neg, depending on forward or backward run
#'diffGMT' can be a vector with the same length as 'runtt'--the DIFFERENCE from GMT at each timestep, which can be dependent on
#	the longitude of current timestep--this is important if need this function to return LOCAL TIME
#Returns MATRIX with following columns:   1) yr 2) mon  3) day  4) hour of day  5)day of week
#	Day of week is value of 0~6, denoting Sun~Sat
#1/25/2000 by JCL
#
#  $Id: weekdayhr.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

runtt<-round(runtt/60)   #change to HOURS
deltad<-sign(runtt)*floor(abs(runtt)/24)  #number of days that need to be adjusted
jstdate<-julian(m=mon,d=day,y=yr)  #compute Julian date of starting time 
jdate<-jstdate+deltad  #vector of resulting Julian dates  

deltahr<-sign(runtt)*(abs(runtt)%%24)  #vector of number of hrs that need to be adjusted
if(length(diffGMT)!=1){
	deltahr<-deltahr+sign(diffGMT)*(abs(diffGMT)%%24)   #if want to function to return LT, adjust for difference from GMT
	sel<-(deltahr<=-24);jdate[sel]<-jdate[sel]-1  #move back one day if deltahr is <=-24
	deltahr[sel]<-deltahr[sel]+24  
}
jhr<-hr+deltahr   #vector of resulting hours
sel<-(jhr<0);jhr[sel]<-jhr[sel]+24  #move back one day
jdate[sel]<-jdate[sel]-1
sel<-jhr>23;jhr[sel]<-jhr[sel]-24  #move forward one day
jdate[sel]<-jdate[sel]+1

datelist<-month.day.year(jdate)
dayofweek<-day.of.week(m=datelist$month,d=datelist$day,y=datelist$year)
result<-cbind(datelist$year,datelist$month,datelist$day,jhr,dayofweek);dimnames(result)<-list(NULL,c("yr","mon","day","hr","weekd"))
result[result[,1]<0,1]<-100+result[,1]
return(result)
}

