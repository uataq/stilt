#!/bin/sh

# usage       :   Run this script "setup.sh". 
#  purposes :  This will set up working directories "Copy0", "Copy1" etc. to run the hymodelc executable within
#                     a main run-time directory "exe", allowing for multiple runs on different nodes.
#                     At least the "Copy0" directory for running STILT is required. At the same place (in "exe") 
#                     a directory "bdyfiles" containing  the files ASCDATA.CFG, LANDUSE.ASC, ROUGLEN.ASC
#                     runhymodelc.bat is created.
#
# tk 05.01.2010 - tkoch@bgc-jena.mpg.de

echo "======================================================================="
echo "This script will set up some subdirectories in the STILT modelling path and some files"
echo " It will create also the working directories \"Copy0\", \"Copy1\" etc. to run the hymodelc executable within"
echo "a main run-time directory (like: \"STILT_Exe\"), allowing for multiple runs on different nodes. At least the"
echo "\"Copy0\" and \"Copy1\" directory is required."
echo "======================================================================="
 
STILT_exe_path=STILT_Exe
Output_path=Output  # for Jena the directory name "Data" is normaly used 
stiltR_path=stiltR

echo "The preferences for the main run-time subdirectory is \"STILT_Exe\""
echo "You want confirm this preferences ? (y/n)"
read confirm

if [ $confirm != "y" ] 
  then
  echo "Enter the subdirectory name which you like instead of \"STILT_Exe\""
  read STILT_exe_path
  echo "Now the main run-time directory allowing for multiple runs on different nodes is:"
  echo $STILT_exe_path
  fi
  
  
echo  "The preferences for the Output subdirectory is: \"Output\""
echo "You want confirm this preferences ? (y/n)"
read confirm

if [ $confirm != "y" ] 
  then
  echo "Enter the subdirectory name which you like instead of \"Output\""
  read Output_path
  echo "Now the output directory is:$Output_path"
  fi 
 
echo "Enter the number of subdirectories (\"Copy1\", (\"Copy2\" etc. you want to create - at least \"1\" "
echo "These directories will be created in $STILT_exe_path subdirectory:"
echo "(please notes that that additional the directory Copy0 will be created for test runs)"
read count_nodes
echo "You want to created in $STILT_exe_path the $count_nodes subdirectories (y/n)"
read confirm

if [ $confirm != "y" ] 
  then
  echo "You have selected \"$confirm\" - so you dont want to continue the script !"
  exit 0
  fi 
  
if  [ $count_nodes -lt 1 ]
  then
  echo "ERROR: You have to select at least \"1\" for creation of one subdirectory" 
  exit 0
  fi
  
 echo "Now enter the  full path where the files LANDUSE.ASC, ROULGEN.ASC, ASCDATA.CFG are"
 echo "normaly this is in the trunk path of your working copy of the STILT inlcuding the corresponding " 
 echo "subdirectories like :  <YourLocalWorkingRepository>/.trunk/stilt_hysplit/bdyfiles"
 read bdyfiles_dir
   
 if [ ! -e  ${bdyfiles_dir}/LANDUSE.ASC ]
  then
  echo "ERROR: The file LANDUSE.ASC is not found in $bdyfiles_dir ! Please start this script from the beginning "
  exit 0
  fi
  
 if [ ! -e  ${bdyfiles_dir}/ROUGLEN.ASC ]
  then
  echo "ERROR: The file  ROUGLEN.ASC is not found in $bdyfiles_dir ! Please start this script from the beginning "
  exit 0
  fi
  
  if [ ! -e ${bdyfiles_dir}/ASCDATA.CFG ]
  then
  echo "ERROR: The file ASCDATA.CFG  is not found in $bdyfiles_dir ! Please start this script from the beginning "
  exit 0
  fi
  
  
 #-----------------make some sub-directories---------
 echo "Now the subdirectories $STILT_exe_path and $Output_path will be created in the STILT modelling directory"
 if [ ! -d $STILT_exe_path ]
  then
  mkdir $STILT_exe_path
  fi 

 if [ ! -d $Output_path ]
  then
  mkdir $Output_path
  fi

 #--------------------------copy some files------------
 echo " "
 echo "Now  the files setStiltparam.r, stilt.bsub.sh, stilt.qsub.csh,  00README.TXT will be created"
 cp ${stiltR_path}/setStiltparam.r .
 cp ${stiltR_path}/stilt.bsub.sh .
 cp ${stiltR_path}/stilt.qsub.csh .
 cp ${stiltR_path}/stilt_simple.sh .
 cp ${stiltR_path}/00README.TXT .

 
 #------------------------------------------------------------  
 echo " "
 echo "Now subdirectories in the $STILT_exe_path will be created"
 cd $STILT_exe_path
 
 i=0
while [  $i -le $count_nodes ]
 do 
   echo " "
   echo $i
   #now create the subdirs
   subdir_name=Copy${i}
   echo "now create subdirectory:$subdir_name"
    if [ ! -d $subdir_name   ]
      then
      mkdir $subdir_name
      fi

   echo "now create runhymodelc.bat in $subdir_name:"
   if [ -e  $subdir_name/runhymodelc.bat ]
     then
     rm $subdir_name/runhymodelc.bat
     fi 
     
   touch $subdir_name/runhymodelc.bat
   echo "cd ${STILT_exe_path}/${subdir_name}" >> $subdir_name/runhymodelc.bat
   echo "hymodelc >>! hymodelc.out" >> $subdir_name/runhymodelc.bat
     
   echo "now copy file LANDUSE.ASC in $subdir_name from $bdyfiles_dir"
   cp  $bdyfiles_dir/LANDUSE.ASC $subdir_name
   
   echo "now copy file ROUGLEN.ASC  in $subdir_name from $bdyfiles_dir"
   cp  $bdyfiles_dir/ROUGLEN.ASC  $subdir_name
      
   echo "now copy file ASCDATA.CFG in $subdir_name from $bdyfiles_dir"
   cp  $bdyfiles_dir/ASCDATA.CFG $subdir_name

  i=$((i+1)) # means: i++
 done
 
echo "============================================================" 
echo  "The  $count_nodes sub-directories   \"Copy0\", \"Copy1\" etc. were created with the"
echo  "corresponding runhymodelc.bat, LANDUSE.CFG ,ROUGLEN. ASCDATA.CFG files"
echo "You have to set  the following variables in the file \"setStiltparam.r\" ":
echo "\"rundir\", \"shlibpath\" and \"metpath\" using any editor"
echo "Script setup.sh is finished  OK"
