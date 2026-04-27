# Download ImmunoNX Example Inputs and Outputs

This document describes how to obtain the input and output data from the ImmunoNX protocol.

## Downloading the HCC1395 FASTQs

The HCC1395 breast cancer cell line data was used as example data for the ImmunoNX protocol. It includes FASTQ files 
for the whole exome and RNAseq from the cancer cell line and data for the matched normal HCC1395 lymphoblastoid line 
([Griffith et al. 2015](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004274)).

### Use `wget` to download TAR files from genomedata.org

``` 
wget https://genomedata.org/hcc1395/fastqs/all/Exome_Norm.tar 
wget https://genomedata.org/hcc1395/fastqs/all/Exome_Tumor.tar
wget https://genomedata.org/hcc1395/fastqs/all/RNAseq_Tumor.tar
tar -xvf Exome_Norm.tar 
tar -xvf Exome_Tumor.tar 
tar -xvf RNAseq_Tumor.tar
```
Note that these files are >100GB combined
