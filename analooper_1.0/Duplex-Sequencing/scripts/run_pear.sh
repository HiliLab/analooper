#!/bin/bash

cd $1

~/pear/bin/pear -m 114 -n 114 -j 2  -f $1.filtered_R1.paired.fq -r $1.filtered_R2.paired.fq -o $1

