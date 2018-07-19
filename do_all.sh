#!/bin/bash

#alltextdays=("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29")

for textday in "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29"
do
 FILES=/group_workspaces/jasmin/geoschem/local_users/lsurl/OMI_ASCII/CH2O/2014/201407/OMI-Aura_L2-OMHCHO_2014m07$textday*".txt"
 for f in $FILES
 do
  g=$(basename "$f")
  bsub -q lotus -o test.out -n 1 -W 2:00 "csh test_scia_in.sh 14 07 "$textday" "$textday" "$g
 done
done