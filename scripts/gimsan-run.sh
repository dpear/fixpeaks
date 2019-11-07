#!/bin/bash

# the line used to call gimsan from /usr/local/bin/gimsan_cmdline/bin

perl gimsan_submit_job.pl --f=~/Projects/fixpeaks/gimsan-input/in-L-not-N.fasta --w=5,6,7,8,10,12 --bg=/data/daniela/Projects/fixpeaks/h99genome/revised_h99_genome.fasta --proc=8 --outdir=~/Projects/fixpeaks/gimsan-output/in-L-not-N/
