#---------------------------------------------------------------------------------------------------
#  $Id: unix.r,v 1.4 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

unix<-function(command, intern = TRUE, ignore.stderr = FALSE){
#redefines unix as system() with intern=T as default
if(nchar(Sys.getenv("COMSPEC"))==0)return(system(command = command, intern = intern, ignore.stderr = ignore.stderr)) #unix/linux
if(nchar(Sys.getenv("COMSPEC"))>0)return(system(command = command, intern = intern)) #PC
}
