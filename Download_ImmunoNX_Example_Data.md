# Download ImmunoNX Example Inputs and Outputs

This document describes how to obtain the input and output data from the ImmunoNX protocol.

## Downloading the HCC1395 Input FASTQ files from SRA

| SRA run | Sample | Output FASTQs |
|---|---|---|
| `SRR37849526` | Matched normal (lymphoblastoid line), normal DNA exome | `Exome_Norm_R1.fastq.gz`, `Exome_Norm_R2.fastq.gz` |
| `SRR37849527` | HCC1395 breast cancer line, tumor DNA exome | `Exome_Tumor_R1.fastq.gz`, `Exome_Tumor_R2.fastq.gz` |
| `SRR37849525` | HCC1395 breast cancer line, tumor RNA-seq | `RNAseq_Tumor_R1.fastq.gz`, `RNAseq_Tumor_R2.fastq.gz` |

The SRA Toolkit (`prefetch` and `fasterq-dump`) is used to retrieve and convert the runs to paired FASTQ files. The steps below run inside a Docker container so no local installation is required. Allow for substantial disk space (>100 GB for the FASTQs, plus temporary space used during conversion).

### Start an interactive session with the SRA Toolkit

The `-v $PWD:$PWD` option mounts your current directory into the container so the downloaded files are written to your working directory, and `-w $PWD` sets it as the working directory inside the container.

```
docker pull ncbi/sra-tools:3.4.1

docker run -it -v $PWD:$PWD -w $PWD ncbi/sra-tools:3.4.1 /bin/bash
```

### Download and convert each run

Within the container, prefetch each run and convert it to paired FASTQ files (`--split-files` writes read 1 and read 2 to separate files; `-e` sets the number of threads):

```
for ACC in SRR37849526 SRR37849527 SRR37849525; do
    prefetch $ACC
    fasterq-dump -e 4 $ACC
done

gzip *.fastq
```

### Rename to the ImmunoNX input filenames

```
mv SRR37849526_1.fastq.gz Exome_Norm_R1.fastq.gz
mv SRR37849526_2.fastq.gz Exome_Norm_R2.fastq.gz
mv SRR37849527_1.fastq.gz Exome_Tumor_R1.fastq.gz
mv SRR37849527_2.fastq.gz Exome_Tumor_R2.fastq.gz
mv SRR37849525_1.fastq.gz RNAseq_Tumor_R1.fastq.gz
mv SRR37849525_2.fastq.gz RNAseq_Tumor_R2.fastq.gz
```

## Downloading the ImmunoNX Output files from Zenodo
All the output files from the ImmunoNX pipeline run have been uploaded to Zenodo at 10.5281/zenodo.17809276. 
The only notable differences between the files on Zenodo and the pipeline outputs are that the RNAseq alignment file output from the pipeline is a BAM file, while the file uploaded in Zenodo is a CRAM file; and the alignment files were renamed to make them more intuitive. 

The pipeline outputs are in a TAR file: hcc1395_pipeline_outs.tar.gz 
Once you have downloaded the TAR file from Zenodo, you can extract the contents using the command:

```
tar -xzf hcc1395_pipeline_outs.tar.gz
```

## Downloading the HCC1395 FASTQs from genomedata

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
