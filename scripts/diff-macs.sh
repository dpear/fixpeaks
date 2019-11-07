#!/bin/bash

t1=/home/daniela/Projects/fixpeaks/bdg3_bed/L1_OD1_3bdg.bed
t2=/home/daniela/Projects/fixpeaks/bdg3_bed/qL1_OD1_3bdg.bed
c1=/home/daniela/Projects/fixpeaks/bdg3_bed/WCE_OD1_3bdg.bed
c2=/home/daniela/Projects/fixpeaks/bdg3_bed/WCE_OD1_3bdg.bed
out=/home/daniela/Projects/fixpeaks/L1-macs-diff
prefix=L1-macs-diff

macs2 bdgdiff --t1 $t1 --t2 $t2 --c1 $c1 --c2 $c2 --outdir $out --o-prefix $prefix
