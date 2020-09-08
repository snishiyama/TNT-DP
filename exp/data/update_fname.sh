#!/bin/bash
cd chronset
for file in *.txt
do
  time=`echo $file | grep -E --only-matching "2020\d+"`
  orifname=`echo $file | grep -E --only-matching "_\d+_.*txt$"`
  cp $file "$time$orifname"

done

