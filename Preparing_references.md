# Preparing References README

There are dozens of input files needed for these pipelines.  These include things like reference genomes, aligner indices, gene annotations, and known variants.   This guide aims to be a comprehensive walkthrough of where to find those files and/or how to create them. Code snippets are included for human, using GRCh38 and Ensembl version 113. 

To create these files yourself, follow the links below, each labelled with key pipelines for which they're needed-

* Alignment
* Germline variant calling
* Somatic variant calling (exome/WGS)
* RNAseq (bulk)

Before running any snippets, make a directory where you'd like this cache to be stored, and set the `BASEDIR` variable to point to it:

```
#BASEDIR=/path/to/annotation_data_grch38_ens113
BASEDIR=/storage1/fs1/mgriffit/Active/griffithlab/gc2596/k.singhal/immuno_manuscript/preparing_references/annotation_data_grch38_ens113
```

## Reference genome 
#### fasta 

We will use the 1000 genomes version of the human GRCh38 build. This reference includes extra decoy and HLA sequences in addition to the alternate haplotypes provided from the GRC consortium. 

```
mkdir -p $BASEDIR/reference_genome
wget -O $BASEDIR/reference_genome/all_sequences.fa https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```


#### fasta index 

From docker image `mgibio/samtools-cwl:1.0.0`

```
/opt/samtools/bin/samtools faidx $BASEDIR/reference_genome/all_sequences.fa
```


#### Sequence dictionary 

From docker image `mgibio/picard-cwl:2.18.1`

```
java -jar /opt/picard/picard.jar CreateSequenceDictionary R=$BASEDIR/reference_genome/all_sequences.fa O=$BASEDIR/reference_genome/all_sequences.dict
```

#### Ensembl Synonyms file

We will download the chromAlias file from UCSC genome browser and restrict it to lines in Ensembl

```
wget -O $BASEDIR/reference_genome/chromAlias.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz
zgrep 'ensembl' chromAlias.txt.gz > chromAlias.ensembl.txt
rm chromAlias.txt.gz
```

## Ensembl/VEP annotation files

#### VEP cache 
Select the most recent [docker image from your desired ensembl version](https://hub.docker.com/r/ensemblorg/ensembl-vep/tags).  Our current release uses `ensemblorg/ensembl-vep:release_113.3`
Inside that container, run the following (mount the appropriate volumes if needed):

```
mkdir $BASEDIR/vep_cache
/usr/bin/perl /opt/vep/src/ensembl-vep/INSTALL.pl --CACHEDIR $BASEDIR/vep_cache --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38
zip -r $BASEDIR/vep_cache.zip $BASEDIR/vep_cache
```

If for some reason the desired cache version and docker image version in use do not match, add `--CACHE_VERSION <version>`

#### Genes/Transcripts GTF 


```
mkdir -p $BASEDIR/rna_seq_annotation
cd $BASEDIR/rna_seq_annotation
wget -O $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.gtf.gz https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
gunzip Homo_sapiens.GRCh38.113.gtf.gz
```

The reference uses "chr" prefixes (chr1, chr2, chr3), while the Ensembl gtf does not (1, 2, 3), so we need to convert the gtf: 

``` 
cd $BASEDIR/rna_seq_annotation/
wget https://raw.githubusercontent.com/chrisamiller/convertEnsemblGTF/master/convertEnsemblGTF.pl
mv $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.gtf $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.gtf.orig
perl convertEnsemblGTF.pl $BASEDIR/reference_genome/all_sequences.dict $BASEDIR/vep_cache/homo_sapiens/113_GRCh38/chr_synonyms.txt $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.gtf.orig >$BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.gtf && rm -f $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.gtf.orig 2> convertEnsemblGTF.stdout
```


#### Transcriptome cDNA reference 

Note that we're not going to convert these coordinates from [1,2,3] to [chr1,chr2,chr3], since Kallisto doesn't include coordinates in it's output.  If pseudobams from kallisto are desired, that conversion would need to be done. We will get both the coding (cdna) and non-coding (ncrna) files and combine them.

```
wget -O $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.cdna.all.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget -O $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.ncrna.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz

cat $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.cdna.all.fa.gz $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.ncrna.fa.gz > $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.cdna.ncrna.fa.gz
```


## Aligner indices 

#### bwa mem index 

From docker image: `mgibio/alignment_helper-cwl:2.1.1`

```
mkdir -p $BASEDIR/aligner_indices/bwamem_0.7.15 
cd $BASEDIR/aligner_indices/bwamem_0.7.15
for i in all_sequences.fa all_sequences.fa.fai all_sequences.dict all_sequences.genome;do 
    ln -s $BASEDIR/reference_genome/$i .
done
bwa-mem2 index all_sequences.fa
```


#### STAR-fusion index 
 In the `trinityctat/starfusion` docker container:

```
mkdir -p $BASEDIR/aligner_indices/star_fusion/temp 
cd $BASEDIR/aligner_indices/star_fusion/temp
/usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl --CPU 10 --genome_fa $BASEDIR/reference_genome/all_sequences.fa --gtf $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.gtf --pfam_db current --dfam_db human --output_dir $BASEDIR/aligner_indices/star_fusion
#based on logfile, this was unsuccessful => issue might be solved by entering a temporary directory but running hte rest of the command the same way. 
# temp directory worked. Now just need to figure out which files to zip. Will use our actual starfusion.zip dir as reference; copying to temp_real_star_fusion_dir and unzipping.

```


## Known Variants

#### Known snps and indels 


dbSNP/indels from the Mills/1000G/etc datasets for use with GATK/Mutect
from the `google/cloud-sdk` image:

```
mkdir $BASEDIR/known_variants
gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf $BASEDIR/known_variants 
TODO BGZIP THIS VCF AND TABIX
gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz $BASEDIR/known_variants
gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi $BASEDIR/known_variants
gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz $BASEDIR/known_variants
gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi $BASEDIR/known_variants

```


#### DoCM cancer variants 

Variants from the [Database of Canonical Mutations](http://docm.info) used for hotspot checking in somatic pipelines

```
wget -O $BASEDIR/known_variants/docm.vcf.gz https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/raw/refs/heads/main/reference_inputs/docm.vcf.gz
wget -O $BASEDIR/known_variants/docm.vcf.gz.tbi https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/raw/refs/heads/main/reference_inputs/docm.vcf.gz.tbi
```


#### Gnomad population frequencies 


This command downloads the V3 of the genome [Gnomad Database](https://gnomad.broadinstitute.org) which has the allele frequencies from 70k WGS samples aligned to GRCh38. This file is ~236gb.

```
#Download gnomad V3 VCFs for each chromosome, concatenate, and filter to AF>0.001.
mkdir $BASEDIR/known_variants/gnomad_temp
# Step 1: Download VCFs
for chr in {1..22} X Y; do
  wget https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v3.1.2/gnomad.genomes.v3.1.2.sites.chr${chr}_trimmed_info.vcf.bgz
done

# Step 2: Concatenate bgzipped VCFs (headers only from first file)
bcftools concat -Oz -o gnomad.genomes.v3.1.2.concatenated.vcf.bgz \
  gnomad.genomes.v3.1.2.sites.chr{1..22}_trimmed_info.vcf.bgz \
  gnomad.genomes.v3.1.2.sites.chrX_trimmed_info.vcf.bgz \
  gnomad.genomes.v3.1.2.sites.chrY_trimmed_info.vcf.bgz

# Step 3: Index the concatenated VCF
bcftools index -t gnomad.genomes.v3.1.2.concatenated.vcf.bgz

# Step 4: Restrict to variants with AF > 0.001
bcftools view -i 'INFO/AF[0]>0.001' -Oz -o gnomad_b38_exome.vcf.gz gnomad.genomes.v3.1.2.concatenated.vcf.bgz
bcftools index -t gnomad_b38_exome.vcf.gz

# Step 5: Copy gnomad VCF to knownvariants and optionally delete original folder
cp $BASEDIR/known_variants/gnomad_temp/gnomad_b38_exome.vcf.gz* $BASEDIR/known_variants/

```

Can create a lite version containing the population allele frequencies from the full vcf using using bcftools. Tested using bcftools version 1.9.

```
# TODO FIGURE OUT WHERE WE WANT TO UPLOAD THESE FILES
bcftools annotate -x ^INFO/AF,INFO/AF_afr,INFO/AF_amr,INFO/AF_asj,INFO/AF_eas,INFO/AF_fin,INFO/AF_nfe,INFO/AF_oth,INFO/AF_sas --threads $THREADS --output-type z -o $BASEDIR/known_variants/gnomad-lite.genomes.r3.0.sites.vcf.gz $BASEDIR/known_variants/gnomad.genomes.r3.0.sites.vcf.gz 
```



#### Common SNPs for somalier concordance 


```
wget -O $BASEDIR/known_variants/somalier_GRCh38.vcf.gz https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
```


## Other Transcriptome files


#### Kallisto transcriptome index 

From docker image `mgibio/rnaseq`

```
kallisto index -i $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.cdna.ncrna.fa.kallisto.idx $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.cdna.ncrna.fa.gz

```


#### Kallisto gene translation table 


From docker image `quay.io/biocontainers/bioconductor-biomart:2.62.0--r44hdfd78af_0`

```
cd $BASEDIR/rna_seq_annotation/
Rscript -e 'library(biomaRt);mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl");t2g <- getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart=mart);t2g <- t2g[order(t2g$ensembl_gene_id), ]; write.table(t2g,"ensembl113.transcriptToGene.tsv",sep="\t",row.names=F,quote=F)'
```


#### Refflat genes for QC 

In docker container: `quay.io/biocontainers/ucsc-gtftogenepred:469--h664eb37_1`

```
gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.gtf /dev/stdout | awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' >$BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.refFlat.txt
```


#### Ribosomal genes for QC 


```
cat $BASEDIR/reference_genome/all_sequences.dict >$BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.ribo_intervals
grep 'gene_biotype "rRNA"' $BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.gtf | awk '$3 == "exon"' | cut -f1,4,5,7,9 | perl -lane '/transcript_id "([^"]+)"/ or die "no transcript_id on $.";print join "\t", (@F[0,1,2,3], $1)' | sort -k1V -k2n -k3n >>$BASEDIR/rna_seq_annotation/Homo_sapiens.GRCh38.113.ribo_intervals
```

#### AGFusion database

```
??? Also need a new docker? https://hub.docker.com/r/mgibio/agfusion
```

## Other

#### Illumina adapters 



```
mkdir $BASEDIR/miscellaneous
cd $BASEDIR/miscellaneous
wget -O illumina_multiplex.fa https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/raw/refs/heads/main/reference_inputs/illumina_multiplex.fa
```


#### QC annotation files 



```
cd $BASEDIR/aligner_indices/biscuit_0.3.8
wget http://zwdzwd.io/BISCUITqc/hg38_QC_assets.zip
unzip hg38_QC_assets.zip
```


#### Capture reagents


```
cd $BASEDIR/miscellaneous
wget -O IDT_xGen_Lockdown_Exome_v1.baits.interval_list https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/raw/refs/heads/main/reference_inputs/IDT_xGen_Lockdown_Exome_v1.baits.interval_list
wget -O IDT_xGen_Lockdown_Exome_v1.targets.interval_list https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/raw/refs/heads/main/reference_inputs/IDT_xGen_Lockdown_Exome_v1.targets.interval_list
```


#### Peptides FASTA file

```
cd $BASEDIR/miscellaneous
wget -O Homo_sapiens.GRCh38.pep.all.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

```