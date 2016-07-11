combine.met<-function(year=99,startmonth=1,endmonth=12, metd="fnl",metlib="/deas/group/stilt/Metdata/",outlib=metlib){
#function to combine arl met data (to avoid hymodelc crashing at beginning of month...)
#combines met files, using cat
#'year': 	2 digit year
#'startmonth':	first month of selected period
#'endmonth':	last month of selected period
#'metd':	type of metfile ("edas" or "fnl")
#'metlib':	where metdata are read from
#'outlib':	where metdata are written to
#2002 by chg
#
#  $Id: combine.met.r,v 1.4 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

#call with: 
#combine.met(03,startmonth=1,endmonth=12,metd="fnl")
#combine.met(03,startmonth=1,endmonth=12,metd="edas")

for(mon in startmonth:endmonth){
for(part in 1:2){
day<-8;if(part==2)day<-24

metfiles<-getmetfile(year,mon,day,12,-15*24,metd=metd,doublefiles=F)

if(metd!="fnl")
outfile<-paste(substring(metfiles[2],1,15),substring(metfiles[1],13,17),sep="")
if(metd=="fnl")
outfile<-paste(substring(metfiles[2],1,10),substring(metfiles[1],8,12),sep="")

print(paste("cat ",metlib,metfiles[2]," ",metlib,metfiles[1]," > ",outlib,outfile,sep=""))
unix(paste("cat ",metlib,metfiles[2]," ",metlib,metfiles[1]," > ",outlib,outfile,sep=""))

}#of for part
}#of for month
}


