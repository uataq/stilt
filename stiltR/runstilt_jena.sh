#!/bin/bash

if [[ "`which bsub`" == "" ]]; then
   echo no queueing command \(bsub\) found
   echo run e.g. on galactica
   exit
fi
if [[ ! -d Runs.done ]]; then
   mkdir Runs.done
fi

rm -f num.txt
for n in `seq 1 10`; do
  while [[ -e num.txt ]]; do
    sleep 1
  done
  echo $n > num.txt
  echo 'bsub stilt part '$n
  echo R CMD BATCH stilt.r stilt_${n}.log | bsub
done

#checkjobs: bjobs, lsload
