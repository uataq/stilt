#Calls Trajecmod
#
#  $Id: stilt.qsub.r,v 1.3 2007-06-29 14:18:29 tnehrkor Exp $
#---------------------------------------------------------------------------------------------------

#get all required R fuctions 
sourcepath<-paste(system("pwd",intern=T),"/",sep="")

if(length(grep("Rsc",sourcepath))==0){ #not the right directory; attach different one
  stop("stilt.r: this is not the right sourcepath")
  sourcepath<-"/project/p1229/jel/Rsc/"
}
print(paste("source directory: ",sourcepath,sep=""))

source(paste(sourcepath,"sourceall.r",sep=""))

totpartarg <- NULL #edit this line to split job up into parts
partarg <- NULL #edit this line to split job up into parts
nodeoffset <- NULL #edit this line to have nummodel <- part+nodeoffset
#Call Trajecmod function, store run info
run.info<-Trajecmod(partarg,totpartarg,nodeoffset)

#save run.info to object with date and time in name
#example: ./Runs.done/setStiltparam.Mar..9.17:40:49.2004.r
savename<-gsub(" ",".",date());savename<-substring(savename,4)
assignr(paste("run.info",savename,sep=""),run.info,"./Runs.done/",printTF=T)

