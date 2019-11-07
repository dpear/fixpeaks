#!/bin/bash

echo "***** Please provide the following fields: ******"
echo "SampleID, Tissue, Factor, Condition, Treatment, Replicate, bamReads, ControlID, bamControl, Peaks, PeakCaller."

read -p 'Samplesheet filename (whole path plz): ' outfile

echo -e "SampleID\tTissue\tFactor\tCondition\tTreatment\tReplicate\tbamReads\tControlID\tbamControl\tPeaks\tPeakCaller" > $outfile

cont='y'

while [ $cont = 'y' ]; do

read -p 'SampleID: ' SampleID
read -p 'Tissue:: ' Tissue
read -p 'Factor: ' Factor
read -p 'Condition (mut or wt): ' Condition
read -p 'Treatment: ' Treatment
read -p 'Replicate (1 or 2): ' Replicate
read -p 'BAM Reads (.bam filepath): ' bamReads
read -p 'ControlID:' ControlID
read -p 'BAM Control (.bam filepath): ' bamControl
read -p 'Peaks (.xls filepath from MACS): ' Peaks

echo -e "$SampleID\tt\tf\t$Condition\tt\t$Replicate\t$bamReads\t$ControlID\t$bamControl\t$Peaks\tbed" >> $outfile

read -p 'Add another sample? (y/n): ' cont
done

cat $outfile
