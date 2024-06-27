# Manual Review README

After the ITB tumor board review all candidates will be evaluated as Accept, Reject, or Review. The candidates marked as Accept and Review will be reviewed in two separate contexts, pVACview and IGV, to verify that they are good neoantigen candidates. 

During the manual review process we primarily will use two files: the Annotated Neoantigen Candidates excel sheet produced from the ITB review and the Peptides 51mer sheet which will will generated. 

To generate the Peptides 51mer excel sheet there is two steps: generating the protein fasta and then generating the manual review files.

## Generate Protein Fasta

First, we will generate an annotated fasta file using a tool which will extract protein sequences surrounding the variant. We will generate one fasta with just the accept and review candidates and another with all proposed candidates. We do this because in somes cases the top candidate may not be the best one (e.g. a different transcript is found to be better during review) so we generate the unfiltered result so that one can consider alternatives.

```bash
gzcat $WORKING_BASE/final_results/annotated.expression.vcf.gz | less
export TUMOR_SAMPLE_ID="hcc1395-tumor-exome"

docker run -it --env WORKING_BASE --env TUMOR_SAMPLE_ID -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/pvactools:4.0.1 /bin/bash

cd $WORKING_BASE/

pvacseq generate_protein_fasta \
      -p $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
      --pass-only --mutant-only -d 150 \
      -s ${TUMOR_SAMPLE_ID} \
      --aggregate-report-evaluation {Accept,Review} \
      --input-tsv $WORKING_BASE/itb-review-files/*.tsv \
      $WORKING_BASE/final_results/annotated.expression.vcf.gz \
      25 \
      $WORKING_BASE/generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa
 
pvacseq generate_protein_fasta \
      -p  $WORKING_BASE/final_results/pVACseq/phase_vcf/phased.vcf.gz \
      --pass-only --mutant-only -d 150 \
      -s ${TUMOR_SAMPLE_ID} \
      $WORKING_BASE/final_results/annotated.expression.vcf.gz \
      25 \
      $WORKING_BASE/generate_protein_fasta/all/annotated_filtered.vcf-pass-51mer.fa

exit
```

Once the peptide fastas have been created we can generate out two final review files using some helper scripts. 

## Generate the Peptide order form

```
export PATIENT_ID="hcc1395"

docker pull griffithlab/neoang_scripts
docker run -it --env WORKING_BASE --env PATIENT_ID -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts /bin/bash

cd $WORKING_BASE
mkdir manual_review

python3 /opt/scripts/generate_reviews_files.py -a itb-review-files/*.xlsx -c generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv -classI final_results/pVACseq/mhc_i/*.all_epitopes.aggregated.tsv -classII final_results/pVACseq/mhc_ii/*.all_epitopes.aggregated.tsv -samp $PATIENT_ID -o manual_review/

python3 /opt/scripts/color_peptides51mer.py -p manual_review/*Peptides_51-mer.xlsx -samp $PATIENT_ID -o manual_review/
```

Open colored_peptides51mer.html and copy the table into an excel spreadsheet. The formatting should remain. 

Maybe some screenshots of these spreadsheets.

## pVACview Review

Binding algorithm support, anchor position
Variants compared to CLE pipeline
IGV review
Somatic variant review (SOP)
Proximal variants review (also look for germline variants, neoantigen candidate may have to be changed to account for it)
RNAseq allele review
Transcript isoform review

## IGV screenshot examples -- different examples from different datasets

![JLF-100-066 – Insertion misalignment - DSTN](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/7a81f397-56e0-410b-857e-14985be0f9eb)


## pvacview screenshots -- different examples from different datasets


Further examples:
https://pvactools.readthedocs.io/en/latest/pvacview/pvacseq_module/pvacseq_vignette.html

