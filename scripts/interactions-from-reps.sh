#!/bin/bash

cat reps-pval-pbg/C_OD1.bed | awk '{if($5!="-"){print $4,"Cqs2",$5}}' > wt-interactions.txt
