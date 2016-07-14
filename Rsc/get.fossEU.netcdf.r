get.fossEU.netcdf<-function(fdate,i,j,ires,jres,spec="co2"){

#function to extract surface fluxes at a specific resolution from flux-maps stored in netcdf files
#fdate: a long vector of dates, each date of format c(YYYY,MM,DD,HH,mm), 
#       the total vector concatenates the individual dates.
#i, j:  particle location indices, referencing to the respective grids at resolution ires, jres
#ires, jres:  resolution at which to extract surface flux averages.
#spec: "CO2" or "CO"
#dpath: path to the netcdf file
#loadpath: path to the fortran shared object (.so file)

if(!any(names(getLoadedDLLs())=="getfossEUnetCDF"))dyn.load(paste(shlibpath,"getfossEUnetCDF.so",sep=""))
#print("log")
#print(paste("get.fossEU.netcdf:",spec))
n<-length(i)
#print(n)
# calling Fortran SUBROUTINE IER_Read (fname, n, date, i, j, ires, jres, f)
#  CHARACTER(255), INTENT(IN)  :: fname                     ! filename
#  INTEGER       , INTENT(IN)  :: n                         ! number of points to get
#  INTEGER       , INTENT(IN)  :: date(5,n),       &        ! YY MM DD hh mm
#                                 i(n), j(n),      &        ! indices of cells
#                                 ires(n), jres(n)          ! averaging parameters
#  REAL          , INTENT(OUT) :: f(n)                      ! averaged flux values

#  length of fname MUST be 255
tmp<-"                                                                                          "
tmp<-paste(tmp,"                                                                                          ")
tmp<-paste(tmp,"                                                                                          ")
#fname <- paste(dpath,"netCDF/CO2.2000a.nc",tmp,sep="")
fname <- paste(emissfile[tolower(spec)],tmp,sep="")
#fname <- ifelse(emissfile=="",paste(emisspath,spec,".2000.nc",tmp,sep="") , paste(emisspath,emissfile,tmp,sep=""))
fname <- substring(fname,1,255)
#print(fname)

if(FALSE){#testing
n     <- 2
fdate  <- as.vector(c( c(2007,1,17,15,9), c(2006,4,30,23,59) ))
i     <- as.vector(c(1,18))
j     <- as.vector(c(39,145))
ires  <- as.vector(c(4,4))
jres  <- as.vector(c(1,2))
}#if(FALSE){#testing

fdate  <- as.vector(fdate)
i     <- as.vector(i)
j     <- as.vector(j)
ires  <- as.vector(ires)
jres  <- as.vector(jres)

f     <- single(length=n)
#res <- .Fortran('ier_read',
#                as.character(fname), as.integer(n), as.integer(fdate),
#                as.integer(i), as.integer(j), as.integer(ires), as.integer(jres), f=f)$f
res <- .Fortran('getfossEUnetCDF',
                as.character(fname), n, as.integer(fdate),
                as.integer(i), as.integer(j), as.integer(ires), as.integer(jres), f)[[8]]
return(res[1:n])
}

