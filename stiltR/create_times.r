#script to create object "Times" w/ starting times for HF
#output is a matrix with 4 columns: 
#-fjul (fractional julian day since 1/1/1960)
#-lat (deg N)
#-lon (deg E)
#-altitude (meters above ground)
#
#  $Id: create_times.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

##############
#path to store receptor information (should be the same as in setStiltparam.r)
path<-"/home/gerbig/Rdat/" 

outname<-"Times.hf" #name for object with receptor information
##############

fjul<-((julian(8,1,2002)*8):(julian(8,31,2002)*8))/8
fjul<-c(fjul,fjul[length(fjul)]+(1:7)/8) #get full last day
fjul<-round(fjul,6)
lat<-42.536
lon<--72.172
agl<-30
assignr(outname,cbind(fjul,lat,lon,agl),path=path,printTF=T)
