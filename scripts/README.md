# fixpeaks/scripts
### Brief description of all scripts in this directory

File | Description
--- | ---
**bed2fasta.sh**                | Wrapper for ```bedtools getfasta``` using the h99 genome
**bed-awk.sh**                  | Formats files with columns: chr, start, end. Too basic of a bash command to exist as its own script, really.
**cluster3-peaks-by-genes.sh**  | Takes sums of bedgraph scores for each gene (that have already been generated for each sample and are in the directory: ```sums-peak-by-gene```, clusters based on sum using kmeans clustering: k=3, and calls genes in the higher-summed two clusters "enriched." Saves results in ```sums3-enriched```.
**cluster3.R**                  | R-script that does the clustering for ```cluster3-peak-by-genes.sh```. Usage: ```Rscript cluster3.R in.bed out.bed```
**cluster-peaks-by-genes.sh**  | Takes sums of bedgraph scores for each gene (that have already been generated for each sample and are in the directory: ```sums-peak-by-gene```, clusters based on sum using kmeans clustering: k=1, and calls genes in the higher-summed cluster "enriched." Saves results in ```sums-enriched```.
**cluster.R**                  | R-script that does the clustering for ```cluster-peak-by-genes.sh```. Usage: ```Rscript cluster.R in.bed out.bed```
**collapse-pval-file.sh**      | Takes in a file with columns: gene, tf, score. Where there is a separate line and score for each tf-gene interaction. Returns a file with columns: gene, tf, score, all-tfs. Where all-tfs is a comma-separated list of all the tf's the gene interacts with. Redundant lines.
**diffbind2homerbed.sh**        | Convert columns from diffbind output to homer input for motif finding and fasta sequence extraction. Mainly for my reference.
**diffbind_analysis.R**         | A smattering of code associated with running the R-package *diffbind* for differential peak binding. This includes code used to read in a file with a list of genes and locations, and annotate a data frame with overlapping genes.
**diff-cutoff-macs-bdg.sh**     | Incomplete code to run the MACS program ```bdgpeakcall``` with a different cutoff for each sample. So far this script only finds the average bedgraph height from each sample in ```ls /data/diana/DKS13_YNBChIPseq/Bedgraphs/*norm_50bp_smooth.bedgraph```.
**diff-macs.sh**                | Wrapper for running ```macs2 bdgdiff```.
**downsample.py**               | Toy python code demonstrating how downsampling from any file could be performed. Courtesy of P Fiziev.
**edgeR-analysis.R**            | All code used in the edgeR analysis. For more refined functionality, see ```run-edgeR.R```
**examine-correlations.R**      | 
find-3-motifs.sh**              |
generate_samplesheet.sh**       |
get-fasta-homer.sh**            |
get-genes.sh**                  |
gimsan-run.sh**                 |
histograms-basic.R**            |
histograms-differences.R**      |
hits.txt**                      |
homer-run.sh**                  |
interactions-from-reps.sh**     |
locate-motif.py**               |
locate-N-motif.sh**             |
macs-bedgraphs.sh**             |
macs-c3-bedgraphs.sh**          |
macs-c4-bedgraphs.sh**          |
macs-od5redo.sh**               |
macs_opt.sh**                   |
macs-redoC.sh**                 |
macs-rerun-k-control.sh**       |
make-network.sh**               |
match-motifs.sh**               |
merge-peak-by-gene.sh** .       |
peakanalysis.R**                |
peak-by-gene-pval.R**           |
peak-by-gene-pval.sh**          |
peak-by-gene.sh**               |
plot-replicates.R**             |
pval-merge.R**      |
pval-merge.sh**      |
resample_bams.sh**      |
resample_rerun_macs.sh**      |
scatter.R**      |
sicer-analysis.sh**      |
sum-peak-by-gene.R**      |
sum-peak-by-gene.sh**      |
tmp.sh**      |
Untitled1.ipynb**      |
Untitled.ipynb**      |
whole-merged-analysis.R**      |
