#!/bin/bash
shopt -s globstar

for num in {2..100}
do
  dir=/home2/danr/marxan_mammals/ph_25cap/run_$num
  run_name=run$num
  echo $dir

  filename=runMarxan_25pc.sh
  file=$dir/$filename

  echo "#$ -N "$run_name"_ph" > $file
  echo "#$ -pe threads 1" >> $file
  echo "#$ -l virtual_free=0.5g,h_vmem=0.6g" >> $file
  echo "#$ -j y" >> $file
  echo "#$ -M dan.rosauer@anu.edu.au" >> $file
  echo "#$ -m bae" >> $file
  echo "hostname" >> $file
  echo >> $file
  echo "cd "$dir >> $file
  echo >> $file
  echo "/usr/local/bin/Marxan input_25pc.dat" >> $file

done
