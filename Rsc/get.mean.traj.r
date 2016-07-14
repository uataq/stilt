get.mean.traj<-function(result,hr=1){
#to get trajectory time series data for storage, average all particles, and keep hr-hourly data (hr=multiple of hour, can be <1)
#sum over all particles and divide by number of particles initially used (nparstilt-1), sum over timesteps within each 2-hour interval
#browser()
#plot(result[,"btime"],result[,"swrad"])

btimeh<-ceiling(result[,"btime"]/hr)*hr

#get rid of initial condition
result<-result[,substring(dimnames(result)[[2]],1,2)!="co"]
#get rid of accumulated flux signals
result<-result[,substring(dimnames(result)[[2]],1,2)!="sv"]

#add all extensiv variables
#colsel<-(substring(dimnames(result)[[2]],1,2)=="sv")+(substring(dimnames(result)[[2]],1,2)=="sa")
colsel<-substring(dimnames(result)[[2]],1,1)=="v"
v<-result[,colsel]#this also gets variable "sampt"
gethrs<-function(vec,hr){return(tapply(vec,hr,sum))}
hv<-apply(v,2,gethrs,btimeh)/(nparstilt) #divide by number of particles initially used

#average over all intensiv variables
gethrm<-function(vec,hr){return(tapply(vec,hr,mean))}
rest<-result[,colsel<1]
hrest<-apply(rest,2,gethrm,btimeh)

resultn<-cbind(hrest,hv)
dimnames(resultn)<-list(NULL,dimnames(resultn)[[2]])
return(resultn)
}
