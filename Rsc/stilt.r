#Calls Trajecmod
#
#  $Id: stilt.r,v 1.7 2010/02/10 19:50:28 trn Exp $
#---------------------------------------------------------------------------------------------------

#
#get all required R fuctions
#
if (! exists("sourcepath")) {
   if (file.exists("stiltR"))
      sourcepath <- paste(getwd(), "/stiltR/", sep="")
   else if (file.exists("Rsc"))
      sourcepath <- paste(getwd(), "/Rsc/", sep="")
   else {
      stop('stilt.r: no stiltR or Rsc directory found.')
      quit(status=1)
   }
}
cat("stilt.r: using sourcepath", sourcepath, "\n")

if (! file.exists(paste(sourcepath, "sourceall.r", sep=""))) stop('stilt.r: "sourcepath" is not a valid source path')
source(paste(sourcepath,"sourceall.r",sep=""))

partinfo <- Sys.getenv(c("STILT_PART", "STILT_TOTPART","STILT_OFFSET"), unset = NA)	# get job partitioning info # new (tk 2011/05/15)
if ( (is.na(partinfo[[3]])) && (!is.na(partinfo[[1]])) && (!is.na(partinfo[[2]] )) ) {
    partinfo[[3]] <-0 # offset is set only optional; means here if the other parts are defined set offset to default value 0 
    }

if (any(is.na(partinfo))) {
   partarg    <- NULL
   totpartarg <- NULL
   nodeoffset <- NULL
} else {
   partarg    <- as.integer(partinfo[[1]])
   totpartarg <- as.integer(partinfo[[2]])
   nodeoffset <- as.integer(partinfo[[3]]) # new (tk 2011/05/17)
}
#Call Trajecmod function, store run info
run.info <- Trajecmod(partarg=partarg, totpartarg=totpartarg, nodeoffset=nodeoffset)

#save run.info to object with date and time in name
#example: ./Runs.done/setStiltparam.Mar..9.17:40:49.2004.r
runs.done.dir <- NULL
if (file.exists('./Runs.done')) runs.done.dir <- './Runs.done/'
if (is.null(runs.done.dir) && file.exists(paste(sourcepath,'Runs.done',sep='')))
  runs.done.dir <- paste(sourcepath,'/Runs.done/',sep='')
if (!is.null(runs.done.dir)) {
  savename <- gsub(" ",".",date())
  savename <- substring(savename,4)
  assignr(paste("run.info",savename,sep=""),run.info,runs.done.dir,printTF=T)
} else {
  cat("stilt.r: Runs.done not found in ./ or sourcepath, not saving run.info\n")
}
