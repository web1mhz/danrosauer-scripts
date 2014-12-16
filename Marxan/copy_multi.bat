#!/bin/bash
shopt -s globstar

for dir in run_{2..100}
do 
  echo $dir
  cp -f /home2/danr/marxan_mammals/sp_25cap_10_x_100/run_1/input_25pc.dat $dir
done
