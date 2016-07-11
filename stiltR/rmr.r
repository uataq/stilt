rmr<-function(xname,path){
#see also assignr() and getr()
#similar to exists(), but 
#removes object at stored location (file path/.Rdataxname)
#3/29/04 by CHG
#
#  $Id: rmr.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------
  system(paste("rm -f ",path,".RData",xname,sep=""))
  return()
}
