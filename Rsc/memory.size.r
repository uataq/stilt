memory.size<-function(){
#instead of memory.size in splus
#for linux only
#2/26/03 chg
#
#  $Id: memory.size.r,v 1.3 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

   user<-system("whoami",intern=T)
   proc<-"R.bin"
   if(length(grep("bgc-jena",system("hostname --long",intern=T)))>0)proc<-"R"
   mems<-system(paste("ps -eo %mem,user,fname | grep ",user," | grep ",proc,sep=""),intern=T)
   return(paste("% memory usage by R: ",mems[1],sep=""))
}
