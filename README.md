# spqn_paper

## Data

### GTEx: 

The original data was downloaded from GTEx portal: https://gtexportal.org/home/datasets

The processed data is available at https://zenodo.org/record/5784626#.YcC7lr3MLxR.

### Drosophila bulk RNA-seq

The original count data was downloaded from: http://bowtie-bio.sourceforge.net/recount/pooled/modencodefly_pooledreps_count_table.txt

The original  metadata was downloaded from: http://bowtie-bio.sourceforge.net/recount/pooled/modencodefly_pooled_phenodata.txt

### scRNA-seq

The original count data was downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45719

### PPI
The original HuRI protein-protein interaction data was downloaded from: http://www.interactome-atlas.org/data/HuRI.tsv


### Regulatome
The original regulatome data was curated using the code from: https://github.com/leekgroup/networks_correction/blob/master/shellscripts/get_true_positive.sh.

The curated data is available at https://zenodo.org/record/5784626#.YcC7lr3MLxR.

### TF gene list
The transcription factor gene list was downloaded from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4825693/bin/NIHMS772895-supplement-Table_S2.xlsx

## Script

### Calculate and adjust correlation matrix
cor_est: remove batch effect, calculate the correlation matrix using GTEx data and adjust the correlation matrix

### Visualize the mean-correlation relationship

signal_ridge: ridgeplot showing the distribution of correlations for genes across different expression levels of GTEx data

corplot: boxplot for the relationship of IQR of correlations and expression levels

slope_min_IQR: scatter plot showing the relationship between the IQR of correlation and the expression level of GTEx data

qqplot: Q-Q plot showing the distribution difference of correlations between different expression levels of GTEx data

scRNA: ridgeplot of scRNA-seq example data

### Visualize the bias in co-expression analysis before and after adjusting the correlation matrix
exp_signal: visualize the weighted distribution of expressions of genes, with the weight decided by the number of estimated co-expression partners, where co-expression signals are decided by setting a threshold on the absolute correlation 

graphical_lasso: visualize the expression bias on co-expression network constricted by graphical lasso

### Examine the expression bias in 'true' co-expression network 
PPI: visualize the weighted distribution of expressions of genes, with the weight decided by the number of interacting proteins of the corresponding protein, according to PPI database 

bias_regulatome: visualize the weighted distribution of expressions of genes, with the weight decided by the number of interacting proteins of the corresponding protein, according to regulatome database 

### Evaluation of SpQN
TP_Partition: compare the true positive of co-expression estimation before and after applying SpQN, conditional on the gene expression level, where co-expression signals are decided by setting a threshold on the absolute correlation, and the gene pairs corresponding to the interacting protein pairs in PPI were used as ground truth

TF_TP: compare the TF-related true positive of co-expression estimation before and after applying SpQN, conditional on the gene expression level, where co-expression signals are decided by setting a threshold on the absolute correlation, and the gene pairs corresponding to the interacting protein pairs in PPI were used as ground truth



### Functions
functions: functions used in the visualization, batch effect removal and correlation matrix adjustment

### Exploratory analysis

scatter_nsampple_IQR: scatter plot showing the change IQR of correlations with the change of the number of PCs removed from the data

TF: visualize the number of predicted co-expressions among TF-related edges

drosophila: explore the mean-correlation relationship of drosophola developmental RNA-seq data

VST: examine the mean-correlation relationship after applying VST on the expression data

combat: examine the mean-correlation relationship after applying combat on the expression data

pool_tissues: examine the mean-correlation relationship for the GTEx data concatinating the expression data of three tissues

