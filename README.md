# spqn_paper

## Data

### GTEx: 

The original data was downloaded from GTEx portal: https://gtexportal.org/home/datasets

The processed data is available at XXX.

### Drosophila bulk RNA-seq

The original count data was downloaded from: http://bowtie-bio.sourceforge.net/recount/pooled/modencodefly_pooledreps_count_table.txt

The original  metadata was downloaded from: http://bowtie-bio.sourceforge.net/recount/pooled/modencodefly_pooled_phenodata.txt

### scRNA-seq

The original count data was downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45719

### PPI
The original HuRI protein-protein interaction data was downloaded from: http://www.interactome-atlas.org/data/HuRI.tsv


### Regulatome
The original regulatome data was curated using the code from https://github.com/leekgroup/networks_correction/blob/master/shellscripts/get_true_positive.sh.

The curated data is available at XXX.

## Script

### Adjust correlation matrix
cor_est

### Visualize the mean-correlation relationship

corplot

signal_ridge: ridgeplot showing the distribution of correlations for genes across different expression levels of GTEx data

slope_min_IQR: scatter plot showing the relationship between the IQR of correlation and the expression level of GTEx data

qqplot: Q-Q plot showing the distribution difference of correlations between different expression levels of GTEx data

scRNA: ridgeplot of scRNA-seq example data

### Visualize the bias in co-expression analysis before and after adjusting the correlation matrix
exp_signal


### Examing the expression bias in 'true' co-expression network 
PPI

###
functions
scatter_nsampple_IQR







TF
