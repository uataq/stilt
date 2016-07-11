# This script is for use with the stilt framework - the function only processes one line of the timesname object instead of the whole object.

# This script will integrate the Wiedenmyer biomass burning inventory into the STILT model for CO, written by Scot Miller
# This script assumes that all of the particle location files will be of named "yearxmonthxdayxhourxlatxlonxagl"

biomassburn<-function(timesname,burnpath="/home/smiller/biomass_burning/",endpath,pathname,nhrs,timesrow) {

# Required Variables
# Timesname: The name of the timesname object.  Unlike regular STILT files, this should already be an R object.
# Burnpath: This is the path to the biomass burning files (These files should be modified as in _script2.txt)
# Endpath:  Where the output file gets written to.
# Pathname: Path to the particle location files and timesname object
# Nhrs: This is the number of hours back in time to calculate the biomass burning influence.
# Note: nhrs should be divisible by 24.
# Biomassburnoutmat: specifies the outputmatrix where the biomassburning values should be dumped
# Timesrow: specifies the row of the timesname object to run.
#
#  $Id: trajecbiomassburn2.r,v 1.3 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------
                       

# Begin to Loop - specify the dates to run and create a matrix to input the data

# Now I have to calculate the appropriate footprints for the times that I'm looking at.	
        
	# The following lines seperate out important input parameters such as the day, month, year, hour of day, etc.
	jdfrac<-timesname[timesrow,1]
        julianday<-trunc(timesname[timesrow,1])   # julianday of this model run from timesname object.
        mmonth<-month.day.year(julianday)[[1]]	  # month of the year
        myear<-month.day.year(julianday)[[3]]     # Year
        modeldays<-month.day.year(julianday)[[2]] # Day of the month
        modelhours<-round((jdfrac-julianday)*24)  # Hour of the day
	agl<-round(timesname[timesrow,4])
        mlat<-round(timesname[timesrow,2],digits=2)
        mlon<- -round(timesname[timesrow,3],digits=2)
        ndays<-round(abs(nhrs)/24,digits=0) # This line converts nhrs into number of days and truncates the number of days to an integer.

	# Set parameters that are dependent on nhrs to be used later in the script.
	# Set nhrs to be used in trajecf foottimes.
		# Note: Trajecfoot must calculate influence for each individual julian day in order for the footprint to be compatible with the
		# temporal resolution of the biomass burning inventory.  Hence, the first time step in foottimes goes back to midnight of
		# of that day and the last time step will tack on the fraction of the julian day required such that the footprint goes back to nhrs.
		if(ndays==1) {
			foottimes<-c(0,modelhours,abs(nhrs)) }
		if(ndays==2) {
			foottimes<-c(0,modelhours,modelhours+24,abs(nhrs)) }
		if(ndays==3) {
			foottimes<-c(0,modelhours,modelhours+24,modelhours+48,abs(nhrs)) }
		if(ndays==4) {
			foottimes<-c(0,modelhours,modelhours+24,modelhours+48,modelhours+72,abs(nhrs)) }
		if(ndays==5){
			foottimes<-c(0,modelhours,modelhours+24,modelhours+48,modelhours+72,modelhours+96,abs(nhrs)) }
		if(ndays==6) {
			foottimes<-c(0,modelhours,modelhours+24,modelhours+48,modelhours+72,modelhours+96,modelhours+120,abs(nhrs)) }
		if(ndays==7) {
			foottimes<-c(0,modelhours,modelhours+24,modelhours+48,modelhours+72,modelhours+96,modelhours+120,modelhours+144,abs(nhrs)) }
		# Check to see if the foottime object exists
			print(paste("Trajecfoot time steps to be used for ",myear,"x",ifelse(nchar(as.character(mmonth))==1,"0",""),mmonth,"x",
				ifelse(nchar(as.character(modeldays))==1,"0",""),modeldays,"x",ifelse(nchar(as.character(modelhours))==1,"0","")
				,modelhours,sep=""))
			print(foottimes)
	# Set trajecf[3] dimnames
		if(ndays==1) {
			trajectimenames<-c(paste(julianday),paste(julianday-1)) }
		if(ndays==2) {
			trajectimenames<-c(paste(julianday),paste(julianday-1),paste(julianday-2)) }
		if(ndays==3) {
			trajectimenames<-c(paste(julianday),paste(julianday-1),paste(julianday-2),paste(julianday-3)) }
		if(ndays==4) {
			trajectimenames<-c(paste(julianday),paste(julianday-1),paste(julianday-2),paste(julianday-3),paste(julianday-4)) }
		if(ndays==5) {
			trajectimenames<-c(paste(julianday),paste(julianday-1),paste(julianday-2),paste(julianday-3),paste(julianday-4),paste(julianday-5)) }
		if(ndays==6) {
			trajectimenames<-c(paste(julianday),paste(julianday-1),paste(julianday-2),paste(julianday-3),paste(julianday-4),paste(julianday-5),paste(julianday-6)) }
		if(ndays==7) {
			trajectimenames<-c(paste(julianday),paste(julianday-1),paste(julianday-2),paste(julianday-3),paste(julianday-4),paste(julianday-5),paste(julianday-6),paste(julianday-7)) }
	
	# The following line calls the appropriate particle location file based on the above parameters.
	ident<-paste(myear,"x",ifelse(nchar(as.character(mmonth))==1,"0",""),mmonth,"x",ifelse(nchar(as.character(modeldays))==1,"0",""),modeldays,"x",
        	ifelse(nchar(as.character(modelhours))==1,"0",""),modelhours,"x",mlat,ifelse(nchar(as.character(mlat))==4,"0","")
                ,"Nx0",mlon,"Wx0",ifelse(nchar(as.character(agl))==3,"0",""),agl,sep="")
	# The following line calls trajecfoot.  It will go back the number of hours specified in nhrs and thus in ndays.
	trajecf<-Trajecfoot(ident=ident,pathname=pathname,foottimes=foottimes)
	print("Trajecfoot done for this receptor point.")

#Pair the biomassburning database down to matching values
        # First pair the biomassburning database only to the required julian days (from present day to number of days that trajecfoot goes backward).
	biomassburningtime<-biomassburning[which(biomassburning[,3]<= julianday & biomassburning[,3] >= (julianday-ndays)),]
        biomassburninglat<-biomassburningtime[which(biomassburningtime[,1] >= 11 & biomassburningtime[,1] < 65),]
	biomassburningfinal<-biomassburninglat[which(biomassburninglat[,2] >= -145 & biomassburninglat[,2] < -51),]

# Now I need to use match to pair up similar points on the matrix.  I'll get 2 vectors and I can multiply them together.
        # First I need to add a julian time name to layers of the array
        	dimnames(trajecf)[[3]]<-trajectimenames
        # Pair the trajecf array down to the matching values
                trajecvector<-vector(length=dim(biomassburningfinal)[[1]])
                for(row in 1:dim(biomassburningfinal)[[1]]) {
                trajecvector[row]<-trajecf[match(round(biomassburningfinal[row,1],digits=6),round(as.numeric(dimnames(trajecf)[[1]]),digits=6)),
                	match(round(biomassburningfinal[row,2],digits=6),round(as.numeric(dimnames(trajecf)[[2]]),digits=6)),
                        match(round(biomassburningfinal[row,3],digits=6),round(as.numeric(dimnames(trajecf)[[3]]),digits=6))]
                } # for row loop
        #Now I can multiply the 2 vectors
        	biomassburningfinal2<-cbind(biomassburningfinal,trajecvector)
                CO1<-biomassburningfinal2[,4]*biomassburningfinal2[,5]*1000
                COppb<-sum(CO1[which(!is.na(CO1))])
		output<-c(timesname[timesrow,1], COppb) # Creates a vector to insert in the biomass burning output matrix.
} # FINISHED
