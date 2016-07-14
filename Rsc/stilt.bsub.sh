#!/bin/bash
#
#---------------------------------------------------------------------------------------------------
#  Purpose
#
#  Script for distributing a STILT run on multiple processors using the LSF queueing system.
#
#---------------------------------------------------------------------------------------------------
#  History
#
#  07/18/2008  Initial version. Stefan Koerner
#
#  $Id: stilt.bsub.sh,v 1.1 2008-08-11 16:08:20 skoerner Exp $
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# example for the addional parameter LF and offset option 
#     ./stilt.bsub.sh  1 1 1 0 Gartow - will run in the  STILT_Exe/Copy1 and produce stilt_01_Gartow.log
#      ./stilt.bsub.sh 1 1 1 1 Helgoland - - will run in the  STILT_Exe/Copy2 and produce stilt_02_Helgoland.log
#-----------------------------------=------------------------------------------

if [[ "$1" == "" || "$2" == "" || "$3" == "" ]]; then
   echo need 3 arguments: start part, end part, total number of parts of the STILT run.
   exit 1
fi

PART1=$1
PART2=$2
TOTPART=$3

#----------now optional define an offset and put additional information to a log-file name------------
LF="" # log file information (additional information in the log file name (e.g. station name
OFFSET=0

COUNT_PARAMETER=$#   # count of caller parameter
if [[ $COUNT_PARAMETER -gt "3" ]]; then
  # add an offset to the stilt run process (for controlling in which STILT_Exe directory the job will be run 
  OFFSET=$4
  fi

if [[ $COUNT_PARAMETER -gt "4" ]]; then
  LF=$5
  LF=_${LF}   # append a underscore
  fi

if [[ "`which bsub`" == "" ]]; then
   echo no queueing command \(bsub\) found
   echo run e.g. on galactica
   exit 1
fi

[[ ! -d Runs.done ]] && mkdir Runs.done

 for P in `seq $PART1 $PART2`; do
   if  [[  ${PART1} == 1 ]]; then
       echo -n submitting STILT part $P of ${TOTPART}"... "
       else
       echo -n submitting STILT part $P of ${TOTPART} starting at ${PART1} "... " 
       fi

   let P_O=$P+$OFFSET
   PF="`printf %2.2i $P_O`"
   bsub -M 2097152 -T 10 -J stilt_${PF} <<- EOD
		export STILT_PART=$P                 # careful: code down to EOD has tab indenting
		export STILT_TOTPART=$TOTPART
                export STILT_OFFSET=$OFFSET  # new and a offset (tk)
		R CMD BATCH --no-restore --no-save --slave stiltR/stilt.r stilt_${PF}${LF}.log
	EOD
done

