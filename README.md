# MERCI: tracing of mitochondrion transfer between single cancer and T cells
MERCI contains two modules:MERCI-mtSNP.py for calling mtSNVs and MERCI R package for predicting the mitochondrial recipient cells and their mitochondrial compositions. 

## MERCI-mtSNP

**MERCI-mtSNP is written in Python3, with the following dependencies:** 

- pandas
- numpy
- pysam
- matplotlib

Please make sure all dependency modules are installed before usiung MERCI-mtSNP.

### Usage
**Type python MERCI-mtSNP.py -h to display all the commandline options:**

| Parameters | Description |
| ------------- | ------------- |
| -h, --help | show this help message and exit |
| -D DATATYPE, --dataType=DATATYPE | The data type of sequencing data. The sequencing data can be of various types, including '10x_scRNA-seq' (default), '10x_mtscATAC-seq', 'smart-seq2', 'bulk_ATAC-seq', 'scATAC-seq', and 'bulk_RNA-seq'.". |
| -o DIRECTORY, --output=DIRECTORY | Output directory for intermediate and final outputs. |
| -S SAMPLEID , --sampleID=SAMPLEID | the sample identifier, also serve as the prefix of output file. if not given, the names of all intermediate or final output files will be automatically set as sampleX. |
| -b PATH_BAM, --Bamfile=PATH_BAM | Input bam file for MT mutation calling |
| -f PATH_FA, --fastafile=PATH_FA | The fasta data of genome reference sequence in the used reference genome file, usually named as XXX.fa |
| -c PATH_BARCODES, --CellBarcode=PATH_BARCODES | This parameter only works for dataTypes with 10x_scRNA-seq or 10x_mtscATAC-seq, the directory where cell barcodes file (barcodes.tsv.gz or barcodes.tsv) generated by cellranger exists |
| -M MQUALITY, --MQcutoff=MQUALITY | The lowest alignment quality that are accepted, the reads with alignment scores below the given value will be discarded, default=5 for scATAC-seq or bulk ATAC-seq, default=255 for other dataTypes |
| -B QCUTOFF, --BQcutoff=QCUTOFF | The base quality cutoff, only alleles with BQ higher than this value will be retained, default=15 for 10x_scRNA-seq, default=25 for other dataTypes |
| -r Species, --ref=Species | This parameter only works for 10x_mtscATAC-seq dataType, user can set 'human' or 'mouse' depending on what species the sequencing data is, default=mouse |
| -l LN, --ln=LN | Only work for 10x_scRNA-seq data type, the maximum length of genomic region for SNP clusters, reads supporting multiple variants within a small genomic region (ln bp) will be removed, default=5 |
| -m MINC, --minC=MINC | For all data types expcept 10x_scRNA-seq, A threshold for coverage, the faction of MT genome that was covered by read counts no less than than this value will be recorded on the generated coverage figure, default=1 |

### Input data format
The processed .bam file generated by alignment software (supporting Cell Ranger, STAR and Bowtie2) is used as the input of MERCI-mtSNP. Make sure the index .bai file is included with .bam file. The fasta file (.fa) in the reference genome folder (essential input for alignment tools) is also needed for sucessfully runing MERCI-mtSNP （controled by -f parameter). For 10x scRNA-seq or mtscATAC0-seq data, the directory of the cell barcodes file generated by Cell Ranger  is  also needed to provide.

Currently MERCI-mtSNP works for 10x scRNA-seq, smart-seq2, RNA-seq, 10x mtscATAC-seq, single cell or bulk ATAC-seq datatypes.

### Example of calling mtSNV from 10x scRNA-seq data:
**Assume the sample id is X, and Cell Ranger alignment output is generated in directory: /X/cellranger/outs. The path of fasta file is: /refdata-cellranger-mm10-2.1.0/fasta/genome.fa, /X/output is the directory where the ouput of MERCI-mtSNP is**
```
python MERCI-mtSNP.py -D 10x_scRNA-seq \
-o /X/output \
-S X \
-b /X/cellranger/outs/possorted_genome_bam.bam \
-f /refdata-cellranger-mm10-2.1.0/fasta/genome.fa \
-c /X/outs/filtered_feature_bc_matrix
```

### Example of calling mtSNV from single cell ATAC-seq data:
**Also assume the sample id is X, and fastq file is aligned using Bowtie2 to the hg38 reference genome. The path of alignment .bam file /X/alignment/X.sort.bam. The path of fasta file is: /UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa**

```
python MERCI-mtSNP.py -D scATAC-seq \ 
-o /X/output \ 
-S X \
-b /X/alignment/X.sort.bam \
-f /UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
```

### Example of calling mtSNV from single cell RNA-seq data with smart-seq2 protocol:
**Assume the sample id is X, and fastq file is aligned using STAR to the hg38 human reference genome. The path of alignment .bam file /X/alignment/X.out.bam. The path of fasta file is: /GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa**

```
python MERCI-mtSNP.py -D smart-seq2 \ 
-o /X/output \ 
-S X \
-b /X/alignment/X.out.bam \
-f /GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

### Example of calling mtSNV from 10x mtscATAC-seq(PMID:32788668) data:
**Assume the sample id is X, and fastq file is aligned using Cell Ranger ATAC pipeline to the mm10 mouse reference genome. The output of Cell Ranger ATAC is generated in directory: /X/cellrangerATAC/outs. 
The path of fasta file is: /refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa**

```
python MERCI-mtSNP.py -D 10x_mtscATAC-seq \ 
-o /X/output \ 
-S X \
-b /X/cellrangerATAC/outs/possorted_bam.bam \
-f /GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
-c /X/cellrangerATAC/outs/filtered_peak_bc_matrix \
-r mouse
```


___*Note: If the folder of .bam file does not contain its .bai index file, the user needs first to generate .bai file. For example using a simple command of samtools to create the bam index file: samtools index *.bam.___

### Output files:
The output directory contains two main output files: *.MT_variants.txt and *.MT_Coverage.csv file (*.Coverage_Cell.csv for 10x scRNA-seq or mtscATAC-seq data). The *.MT_variants.txt contains the annotated information of retrieved mtSNVs, *.MT_Coverage.csv or *.Coverage_Cell.csv file records the coverage information in mitochondrial genome for each cell or sample.


## MERCI R package
### Install MERCI package
```
library(devtools)  
install_github(repo='shyhihihi/MERCI/MERCI')  
```
___*Note: When the system reminds you to update dependent packages, we recommend not to update___
### Usage
**This text will use a benchmark example data to illustrate how to use MERCI package**  
After the MERCI package is successfully installed, load the package.

`library(MERCI)`

Read the file of mitochondrial variants (*.MT_variants.txt’) called from MERCI-mtSNP into R environment. Here I prepared an example data ‘example.MT_variants.txt’. (Users can find the dataset of example files in https://github.com/shyhihihi/MERCI/tree/main/data)
```
varFile  <- './example.MT_variants.txt' ;  
mtSNV_table <- readMTvar_10x(varFile, minReads=1000) ;  
```
Here, varFile is the path of '*.MT_variants.txt' file, and minReads=1000 indicates that only Cell with MT reads < 1000 will be filtered out. readMTvar_10x  function works only for the variant file of 10x scRNA-seq data. If the user used other datatypes, such as ATAC-seq, smart-seq2 or bulk RNA-seq, etc., please use the code below:
```
varFile  <- '. / XXX.MT_variants.txt'  
MT_variants <- readMTvar(varFile, cellname = "XXX") 
```

### MERCI LOO for real-world application
**We first show how to run MERCI LOO pipeline to predict mitochondrial receiver cells, which means there is no available reference data of donor and pure non-receiver cells. If the user have reference data, we recommend to use MERCI regular pipeline (see the last section)**  
Assuming T cells are mitochondrial donor cells and the cancer cell population is a mixture of receivers and non-receivers. In this example data, we mixed 300 CC (receiver) and 500 MC (non-receiver) cancer cells together. Load the cell information data:
```
load('. /cell_info.RData')  
T_cells <- cell_info$cell_name[cell_info$cell_type=='T cell']  
Cancer_cells <- cell_info$cell_name[cell_info$cell_type=='cancer cell']  
```

Read the MT read-coverage file generated by MERCI-mtSNP into R.
```
CoverFile <- './example.Coverage_Cell.csv'  
selected_Cells <- c(T_cells, Cancer_cells)  
MTcoverage_inf <- readCoverage_10x(CoverFile, S.cells=selected_Cells)  
```

Also, readCoverage_10x function works only for the coverage file of 10x scRNA-seq data. If the user used other datatypes, please use the code below to read the coverage information:
```
CoverFile <- './XXX.MT_Coverage.csv'  
MT_cov <- readCoverage(CoverFile)  
```

Preprocess the MT mutation data and generate a mtSNV vaf matrix (variants * cells)
```
s.mtSNV_table <- mtSNV_table[mtSNV_table$Cell%in%selected_Cells, ]  
mtSNV_ma <- MTmutMatrix_transform(MT_variants=s.mtSNV_table, MT_coverage=MTcoverage_inf, donor_cells=T_cells, mixed_cells=Cancer_cells, min_d=5, min_observeRate= 0.1)   
```

We focus on the T cells (potential donor cells) and cancer cells (potential receiver cells).   
For these selected cells, we get the VAF matrix of their mtSNVs (mtSNV_ma).    
Next, we identify the donor cell (T cell) enriched mtSNVs

`s.muts <- Denrich_mtMut_extr(varMatrix=mtSNV_ma, donor_cells=T_cells, mixed_cells=Cancer_cells, OR=2, Nmut_min=2)`

Here, to get features (MT variants) as many as possible for downstream prediction, we only set the Odds Ratio (T cell vs Cancer cell) > 2 and leave  p and q values as defaut alone.  
Next, we calculate the ‘effective count statistic’(Neff)  of donor cell enriched mtSNVs and the DNA ranks for all candidate receivers (mixed cancer cells).  

`DNA_rank <- Cell_Neff_cal(varMatrix=mtSNV_ma, MT_variants=s.mtSNV_table, MT_coverage=MTcoverage_inf, donor_cells=T_cells, mixed_cells=Cancer_cells, mutFeatures=s.muts, adjust=FALSE)`

The parameter ‘adjust’ of Cell_Neff_cal function indicates whether adjusting the count statistics with Neff statistics. when adjust=FALSE, then the original count of observed donor-enriched mtSNV in each candidate receiver cell will be used to calculate the DNA rank. Otherwise, Neff will be used as the count statistic for DNA rank calculation, default=TRUE. Notably, When using Neff to calculate DNA rank score, it will take a long time to get the results depending on the number of ‘mixed_cells’ and the number of ‘mutFeatures’.  

Next using MERCI LOO pipeline to estiamte the relative abundance of transferred mitochondria (decovolution analysis)and RNA rank for each input cancer cell.

First Load gene expression data (the expression matrix must include the expression porfiles for MT genes),
```
load('./cell_exp.RData')  
library(Matrix)  
cell_exp <- cell_exp[, selected_Cells]  
RNA_rank <- MERCI_LOO_MT_est(cell_exp, reciever_cells=Cancer_cells, donor_cells=T_cells, organism='mouse')  
```
MERCI_LOO_MT_est will return the estimated MT constitutes for all candidate receiver cells and the transformed RNA ranks. The parameter ‘organism’ must be assigned accurately, currently, we only support Human and mouse species.  
The output of MERCI_LOO_MT_est function (here the R object 'RNA_rank') contains the estimated relative fraction of T cell-derived mitochondria (the first column ' Donor_MT_ind') for each input cancer cell.  
It is easy to find that the relative faction of T cell-derived MT in the real receiver cells (CC cancer cells) is much higher than non-receiver cells (MC cancer cells), see R codes below.
```
ture_receiver_cells <- cell_info$cell_name[cell_info$culture_history=='CC_cancer'] ;
non_receiver_cells <-  cell_info$cell_name[cell_info$culture_history=='MC_cancer'] ;
boxplot(RNA_rank[ture_receiver_cells, 'Donor_MT_ind'], RNA_rank[non_receiver_cells, 'Donor_MT_ind'], main='The relative abundance of foreign MT', names=c('CC', 'MC') ) 
```
![Image text]( https://github.com/shyhihihi/MERCI/blob/main/images/boxplot_CCvsMC.jpeg)


Significance estimation to test if true-receivers are included in the input cancer cell population based on Rcm values.

`CellN_stat <- CellNumber_test(DNA_rank, RNA_rank, Number_R=1000)`

The statistic Rcm value will be returned. If there is Rcm >1 at any cutoff, this means receivers are high likely to be sufficiently included in the input mixed cells. Let’s look at the results for the example data:

![Image text]( https://github.com/shyhihihi/MERCI/blob/main/images/MERCIv2_LOO_Rcm.jpeg)

For rank cutoffs at top rank 10-80%, the Rcm values are consistently > 1. This means the captured number of positive calls is significant (p < 0.001) and non-random, true receivers are thus considered sufficient in the input cancer cells.  
Let's explore the results if there is no ture-receiver within the input cancer cells.
```
test.CellN_stat <- CellNumber_test(DNA_rank[non_receiver_cells, ], RNA_rank[non_receiver_cells, ], Number_R=1000)
```
![Image text]( https://github.com/shyhihihi/MERCI/blob/main/images/Rcm_nonReceivers.jpeg)

As the figure above shows, there is no Rcm value >1 at any rank cutoff. This means it is highly possible that no real receivers (or no sufficient receivers) are included in the input cell set.

We next selected a rank cutoff to predict the MT receivers. Here, we used the cutoff at top rank 50%, which is a good choice to balance the sensitivity, specificity and precision.

`MTreceiver_pre <- MERCI_ReceiverPre(DNA_rank, RNA_rank, top_rank=50)`

**The prediction task is completed here.**  
Let’s look at the performance of prediction results.
```
t.stat <- table(cell_info[Cancer_cells, 'culture_history'], MTreceiver_pre[Cancer_cells, 'prediction'])  
t.stat
```
![Image text](https://github.com/shyhihihi/MERCI/blob/main/images/MERCIv2_t.stat.jpg)

According to the equations blelow:  
Precision=TP/(TP+FP)=226/(226+41)=84.6%  
Specificity=TN/(TN+FP)=459/(459+41)=91.8%  
Sensitivity=TP/(TP+FN)=226/(226+74)=75.3%  
The precision (also called positive predictive value) reached > 84%. And the specificity and  the sensitivity are 92% and 75%, respectively. If we selected a more rigorous cutoff (e.g. top 40% or higher), the precision and specificity will increase at the cost of reduced sensitivity.

### MERCI regular
**If there is reference data provided, the prediction performance will be better when using MERCI regular pipeline conpared to MERCI LOO pipeline. Thus we recommond to use regular MERCI pipeline as isllustrated below:**  
Load the independent reference data for pure non-receivers of cancer cells (additional MC cells), including the cell annotation data and gene expression data.
```
load(‘./cell_info_nonReceivers.RData’)  
load(‘./cell_exp_nonReceivers.RData’)  
ref_noRe_cells <- cell_info_nonReceivers$cell_name  
selected_Cells <- c(T_cells, Cancer_cells, ref_noRe_cells)  
c.genes <- intersect(rownames(cell_exp), rownames(cell_exp_nonReceivers))
cell_exp2 <- cbind(cell_exp[c.genes, ], cell_exp_nonReceivers[c.genes, ])
```
Read the file of mitochondrial read-coverage into R environment based on selected cells, and generate the corresponding vaf matrix for mtSNVs.
```
MTcoverage_inf <- readCoverage_10x(CoverFile, S.cells=selected_Cells)  
s.mtSNV_table <- mtSNV_table[mtSNV_table$Cell%in%selected_Cells, ]  
mtSNV_ma2 <- MTmutMatrix_transform(MT_variants=s.mtSNV_table, MT_coverage=MTcoverage_inf, donor_cells=T_cells, mixed_cells=Cancer_cells, Ref_nonReceivers = ref_noRe_cells, min_d=5, min_observeRate= 0.1)  
```

Calculated the DNA and RNA ranks for the input mixed cells of cancer (mixed_cells), using the data of T cells and new loaded non-receivers as reference.
```
s.muts <- Denrich_mtMut_extr(varMatrix=mtSNV_ma2, donor_cells=T_cells, mixed_cells=Cancer_cells, Ref_nonReceivers=ref_noRe_cells, OR=2, Nmut_min=2)  
DNA_rank2 <- Cell_Neff_cal(varMatrix=mtSNV_ma2, MT_variants=s.mtSNV_table, MT_coverage=MTcoverage_inf, donor_cells=T_cells, mixed_cells=Cancer_cells, Ref_nonReceivers = ref_noRe_cells, mutFeatures=s.muts, adjust=FALSE)  
RNA_rank2 <- MERCI_MT_est(cell_exp2, mixed_cells=Cancer_cells, donor_cells=T_cells, Ref_nonReceivers=ref_noRe_cells, organism='mouse')  
```

Also, perform significance estimation first to obtain the Rcm statistics.

`CellN_stat2 <- CellNumber_test(DNA_rank2, RNA_rank2, Number_R=1000)`
![Image text]( https://github.com/shyhihihi/MERCI/blob/main/images/MERCIv2_regular_Rcm.jpeg)  

Rcm values are consistent > 1 at cutoffs from top rank 10-80%. We next used the same cutoff 50% to predict the mitochondrial receivers.
```
MTreceiver_pre2 <- MERCI_ReceiverPre(DNA_rank2, RNA_rank2, top_rank=50)  
t.stat2 <- table(cell_info[Cancer_cells, 'culture_history'], MTreceiver_pre2[Cancer_cells, 'prediction'])
t.stat2  
```
![Image text](https://github.com/shyhihihi/MERCI/blob/main/images/MERCIv2_t.stat2.jpg)


Compared to the results of prediction without reference data (the results of MERCI LOO pipeline), we can easily find the performance improved significantly with sensitivity = 79%. At the same time, there is nearly no change of precision (85%) and specificity (91%). It is ok for using the MERCI LOO pipeline if the user does not have additional reference data, but MERCI regular are highly recommended if there is reference data for pure non-receivers and pure donor cells.

