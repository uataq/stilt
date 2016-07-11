assignr<-function(xname, value, path="", printTF=FALSE,gz=F){
#see also getr() and existsr()
#to avoid large databases that take too long to load memory (as R tries...)
#similar to assign(), but object is removed after call to assignr
#assigns object to name xname, saves it with save() under path/.RDataxname
#and REMOVES it from local database (search()[1])
#example: assignr("test",temp,path="/mydir/") creates file "/mydir/.RDatatest"
#that can be attached and contains a single object with name "test"
#printTF can be set to TRUE to print a statement that the object has been assigned
#2/27/04 by CHG
#modified to allow compression of files (flag gz)
#
#  $Id: assignr.r,v 1.4 2008-08-12 08:51:22 skoerner Exp $
#---------------------------------------------------------------------------------------------------
     assign(x=xname,value=value,pos=1)
     save(list=xname,file=paste(path,".RData",xname,sep=""))
     remove(list=xname,pos=1)
     if(gz){
       unix(paste("gzip ",path,".RData",xname,sep=""))
       xname<-paste(xname,".gz",sep="")
     }
     if (printTF) cat(path, ".RData", xname, " created\n", sep="")
}
