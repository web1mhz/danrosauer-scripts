#!/bin/bash
shopt -s globstar

for num in {1..100}
do
  dir=/home2/danr/marxan_mammals/sp_25cap_10pc/run_$num

  filename=runMarxan_25pc.sh
  file=$dir/$filename
  echo $file

  qsub $file

done
