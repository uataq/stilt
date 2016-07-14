#!/bin/csh -f

set R_cmd=`which R`
if ($R_cmd[1] == "R:") then
  set R_cmd=$R_HOME/bin/R
  echo "R not in path, using $R_HOME to define R_cmd=$R_cmd"
endif

unset totpart
unset parts
set setup=0
set dryrun=0
set workdir=`pwd`
set hostlocal=`hostname`
set wall="168:00:00"
set exename='Exe'
if ($?STILT_HOME) then
 set stiltdir=$STILT_HOME
else
 set stiltdir=`pwd`
endif

unset qsub
switch ($hostlocal)
case humboldt*:
  echo "Using default settings for humboldt"
  set qsub_def="$HOME/bin/qsub_hum -lwalltime=$wall"
  breaksw
case beehive*:
  echo "Using default settings for beehive"
  set qsub_def="qsub -pe mpich 1 $HOME/bin/qsub_generic.csh"
  breaksw
default:
  echo "No match found for host=$hostlocal, using generic default settings"
  set qsub_def="qsub"
endsw

# Parse arguments for options:
set dashdone=0
while ($#argv > 0 && $dashdone == 0) 
   switch ($argv[1])
   case -help:
    goto usage
    breaksw
   case -s:
    set setup=1
    shift argv
    set dashdone=0
    breaksw
   case -n:
    set dryrun=1
    shift argv
    set dashdone=0
    breaksw
   case -q:
    shift argv
    set qsub="$argv[1]"
    shift argv
    set dashdone=0
    breaksw
   case -w:
    shift argv
    set wall="$argv[1]"
    shift argv
    set dashdone=0
    breaksw
   case -o:
    shift argv
    set nodeoffset="$argv[1]"
    shift argv
    set dashdone=0
    breaksw
   case -e:
    shift argv
    set exename="$argv[1]"
    shift argv
    set dashdone=0
    breaksw
   case -t:
    shift argv
    set totpart=$argv[1]
    shift argv
    set dashdone=0
    breaksw
   case -p:
    shift argv
    if ($?parts == 0) then
       set parts=$argv[1]
    else
       set parts=( $parts $argv[1] )
    endif
    shift argv
    set dashdone=0
    breaksw
   default:
    set dashdone=1
   endsw
end

#repeat here to let $wall take effect for qsub_def:
switch ($hostlocal)
case humboldt*:
  echo "Using default settings for humboldt"
  set qsub_def="$HOME/bin/qsub_hum -lwalltime=$wall"
  breaksw
case beehive*:
  echo "Using default settings for beehive"
  set qsub_def="qsub -pe mpich 1 $HOME/bin/qsub_generic.csh"
  breaksw
default:
  echo "No match found for host=$hostlocal, using generic default settings"
  set qsub_def="qsub"
endsw

if ($?qsub == 0) then
   set qsub="$qsub_def"
endif

if ($#argv > 0) then
   echo "$0 : unrecognized arguments: $*"
   goto usage
endif

if ($?totpart == 0) goto usage
if ($?parts == 0) then
   @ num=0
   while ( $num < $totpart )
     @ num++
    if ($?parts == 0) then
       set parts=$num
    else
       set parts=( $parts $num )
    endif
  end
endif
if ($?nodeoffset == 0) then
   set nodeoffset=NULL
endif
if ($setup) then
   if (-d $stiltdir/Exe/Copy1 == 0) then
      echo "Cannot find $stiltdir/Exe/Copy1, rerun with proper setting for STILT_HOME"
      exit 1
   endif
   set lnfiles=( `find $stiltdir/${exename}/Copy1 -type l -print` )
endif

foreach num ( $parts )
  if ($nodeoffset == "NULL" || $nodeoffset == 0) then
     set exenum=$num
  else
     @ exenum=$num
     @ exenum+=$nodeoffset
  endif
  if ($setup) then
   set olddir=`pwd`
   if (-d ${exename}/Copy$exenum == 0) mkdir -p ${exename}/Copy$exenum
    cd ${exename}/Copy$exenum
    echo -n "Checking link files in:";pwd
    foreach tst ( $lnfiles )
      if (-e $tst:t == 0 || -z $tst:t) then
         echo "Doing: ln -sf $tst ."
         if ($dryrun == 0) ln -sf $tst .
      endif
    end
    cd $olddir
  endif
  set stiltr_dir=${stiltdir}/stiltR
  if (-d $stiltr_dir == 0) set stiltr_dir=${stiltdir}/Rsc
  set stiltr_in=${stiltr_dir}/stilt.r
  if (-e $stiltr_in == 0) then
     echo "Cannot find stilt.r in either Rsc or stiltR subdirectories of $stiltdir"
     echo "Rerun with proper setting of STILT_HOME"
     exit 1
  endif 
  if (-e stiltR == 0) then
     echo "Doing: ln -sf $stiltr_dir stiltR"
     if ($dryrun == 0) ln -sf $stiltr_dir stiltR
  endif
  sed -e "/^ *totpartarg/s/NULL/$totpart/" -e "/^ *partarg/s/NULL/$num/" \
    -e "/^ *nodeoffset/s/NULL/$nodeoffset/" ${stiltr_in} >! stilt.r.$exenum
  echo "Submitting $num of $totpart"
  echo '#\!/bin/csh -f' >! stiltr.job.$exenum
  echo "$R_cmd --no-save < stilt.r.$exenum >&! stilt.r.$exenum.out" >> stiltr.job.$exenum 
  chmod +x stiltr.job.$exenum
  echo "$qsub $workdir/stiltr.job.$exenum"
  if ($dryrun == 0) $qsub $workdir/stiltr.job.$exenum
end

exit

usage:
  echo "Help message for $0"
  echo "This script submits multiple stiltr jobs via qsub by editing part,totpart,nodeoffset in stilt.r"
  echo "Usage: $0 options"
  echo "Options:"
  echo "Required: -t <totpart>: total number of parts to break the job into"
  echo "Optional: -p <part> [-p ...] : part(s) to run; default: run all (1-totpart)"
  echo "Optional: -o <nodeoffset>: number to add to part to determine ${exename}/Copyn directory (default:0)"
  echo "Optional: -q <qsub>: qsub command; default: $qsub_def"
  echo "Optional: -w <wall>: walltime string (hh:mm:ss) to be used for default humboldt qsub command; default: $wall"
  echo "Optional: -s : setup Copy directories if needed: replicate symbolic links from Copy1"
  echo "Optional: -n : dry-run only, do not issue commands for setup or qsub"
  exit 1
