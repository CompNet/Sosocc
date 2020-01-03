#!/bin/bash

supplied=1
if [ $# -eq 0 ]; then
  echo "No arguments supplied"
  supplied=0
fi

OUTFILENAME="used-ram-memory.txt"
if [ $supplied -eq 1 ]; then
  OUTFILENAME="$1"/"$OUTFILENAME"
fi

echo 'output filename is ' $OUTFILENAME

# ----------------------------------


# we wait 5 secs (chosen arbitrarly), becase we need that the partitioning process has already started to get its pid
# since this scripts are executed before that, it should work
sleep 5;

# retreive the parent pid:
# first, ensure that the process is executed from R
pProcString=`pstree -lps -N user | grep R | grep sh | grep java`
# this gives something line: "bash(23956)---R(18843)---sh(29442)---java(29443)-+-{java}(29444)"
# for instance, the parent pid is 29443, so retreive it

delim='---'
arrIN=(${pProcString//$delim/ })
delim2='-+-'
arrIN2=(${arrIN[3]//$delim2/ })
jppid=`echo ${arrIN2[0]} | sed 's/.$//' | sed 's/^java(//'`
echo $jppid

# -----------------------------------------------------------------


while true; do

  jpstatus=`cat /proc/$jppid/smaps`
  # you can inspect $?, which contains the exit code of the last command (zero for success, non-zero for failure)
  if [[ $? != 0 ]]; then
	break # exit the infinite loop
	echo 'no partitioning process found. Exiting ..'
  fi

  # get memory consumption in GB
  # newVal=`cat /proc/$jppid/smaps | awk '/Pss/ {mem += $2} END {print mem, "kB"}' | awk '{ byte =$1 /1024/1024; print byte }'` # in gB
  newVal=`cat /proc/$jppid/smaps | awk '/Pss/ {mem += $2} END {print mem}'` # in kB
  if [ ! -f $OUTFILENAME ]; then
    echo $newVal > $OUTFILENAME
  else
    currVal=`cat $OUTFILENAME`
    if [ $newVal -gt $currVal ]; then
      echo $newVal > $OUTFILENAME
    fi
  fi

  sleep 1;
done

# do it in R
# pid = system("pgrep record-mem.sh", intern=TRUE)
# system(paste("kill -9", pid)) 
