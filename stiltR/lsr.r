lsr<-function(path,pattern=""){
#lists all r objects in 'path'
#
#  $Id: lsr.r,v 1.4 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

tmp<-dir(path,pattern=paste(".RData",pattern,sep=""),all.files=TRUE)
return(substring(tmp,nchar(paste(".RData",sep=""))+1))
}
