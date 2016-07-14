read.bg<-function(spec,datename,pathin,pathout){
#reads initial field for CO2, CO, CH4, or H2, created using sas

#spec is any subset of c("co","co2","ch4","h2")
#datename is a string indicating time range, e.g. "1_1_99_12_31_02" for 1/1/1999 to 12/31/2002
#pathin is path where to read ASCII data from
#pathout is path where to assign boundary condition objects

#OUTPUT
#assigns *.ini objects (3-D array altitude * Latitude * Day since 1/1/1960)
# e.g. read.bg(spec=c("CO2","CO"),datename="1_1_99_12_31_02",pathin="/group/stilt/Boundary/",pathout="~johnlin/Rdat/")
#
#  $Id: read.bground.r,v 1.3 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

#sasdate is 0 at 1/1/1960 (i.e. elapsed days since 1/1/1960)
for(sname in tolower(spec)){
bg<-read.asc(paste(pathin,sname,"bg",datename,".asc",sep=""))
print(dim(bg))
bg<-bg[order(bg[,"sasdate"],bg[,"lat"],bg[,"agl"]),]
print(dim(bg))

bg.ini<-array(bg[,paste(sname,"m",sep="")],dim=c(length(unique(bg[,"agl"])),length(unique(bg[,"lat"])),length(unique(bg[,"sasdate"]))))
print(dim(bg.ini))
dimnames(bg.ini)<-list(as.character(unique(bg[,"agl"])),as.character(unique(bg[,"lat"])),as.character(unique(bg[,"sasdate"])))
assignr(paste(sname,".ini",sep=""),bg.ini,path=pathout,printTF=T)
mdy1<-month.day.year(min(bg[,"sasdate"]))
mdy2<-month.day.year(max(bg[,"sasdate"]))

print(paste(sname," valid from ",mdy1[1],"/",mdy1[2],"/",mdy1[3]," to ",mdy2[1],"/",mdy2[2],"/",mdy2[3]," (m/d/y)",sep=""))
}
}
