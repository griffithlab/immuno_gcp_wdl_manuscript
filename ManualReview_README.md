# Manual Review README

The Immunotherapy Tumor Board (ITB) reviews all candidates and evaluates them as Accept, Reject, or Review. The candidates marked as Accept and Review will be reviewed in this 'Manual Review' section in two separate contexts, pVACview and IGV, to verify that they are good neoantigen candidates. 

During the manual review process, we will primarily use two files: the Annotated Neoantigen Candidates excel sheet produced from the ITB review and the Peptides 51mer sheet which we will generate below. 
The Peptides 51mer sheet is generated in a certain format specified by the peptide manufacturer used for the neoantigen vaccine trials.

To generate the Peptides 51mer excel sheet there are two steps: (1) generating the protein fasta and (2) generating the manual review files.

## Generate Protein Fasta

First, we will generate an annotated fasta file using a tool that will extract protein sequences surrounding the variant. We will generate one fasta file with just the 'accept' and 'review' candidates and another with all proposed candidates. We do this because in some cases, the top candidate may not be the best one (e.g. a different transcript is found to be better during manual review) so we generate the unfiltered result so that one can consider alternatives.

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

Once the peptide fastas have been created we can generate two final review files using some helper scripts. 

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

## Reviewing QC data
- running the qc scripts?
      - the qc cutoffs
- end bias examples
![5120-30 end-bias horseshoe GOOD]https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/2cea9675-5f30-45a2-b6ad-8afe458c7099)

## Reviewing HLA alleles

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

### Unaccounted for Germline Variants
(Leidos Case 5120-28) MST1R (screenshots below) had an upstream frameshift mutation which was not accounted for in the 51mer sequences. The sequence has been manually fixed and the germline variant is marked in purple. We have not yet determined why this germline variant was not automatically accounted for. 

Mutation on the left exon and germline variant on the right exon.
![5120-28 MST1R germline variant](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/59aeae0f-583a-4ad9-958c-af6aeb20c35c)

A -> C heterozygous germline variant results in Gln (Q - CAA) becoming  Pro (P - CCA)
![5120-28 MST1R germline variant AA sequence](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/9390daaa-9890-4863-821f-b44791c001cc)

The variant happens in phase with the germline variant only some of the times.


A really detailed description on exploring a germline varaint in 5120-18

![5120-28 MST1R germline variant no in phase with variant](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/5f2075b2-6ffc-44b2-9efc-b3ad1fcf5ca6)

### Multiple RNA Splicing
(Leidos Case 5120-27) Just a note that candidate MT.60.CPEB2.ENST00000507071.6.missense.251P/L, has multiple RNA splicing, but all transcripts give rise to the same peptide
![5120-27 CPEB2 alternative splicing](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/d341bc01-f80c-4817-881f-1695380f4085)

![5120-27 CPEB2 alternative splicing peptide sequences](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/a13601c7-10e6-4548-8ec8-50cd00f0d047)

### Framshift miscount

Leidos 5120-19 - RNA counts were corrected for two single base insertions showing that their expression was indeed supported despite initially having RNA VAFs of 0%: HES1 FS46 and RREB1	FS781-782

### Dunicleotide Variant

5120-16: AASDH (V140L): there is a true dinucleotide variant that is captured separately and therefore this one should be rejected since it only addresses one incorrectly.

### Complex Variant
(TODO: need to check if we can include this as an example)
MCDB044: MAP2: There is a DNP adjacent to and in phase with a 28-base deletion. Manually looked at the wild-type sequence, and all mutations and translated the sequence to predict the appropriate peptide sequence that should be incorporated in the design.  

### transcript not expressed
5120-17: ADPRHL1 - was removed from consideration - the expressed transcript does not contain the mutation

### Reference Match KRAS G12R candidate

5120-06
KRAS G12R candidate
According to CancerHotspots, the R amino change is the 4th most common (D > V > C > R).  

Other than the reference match issue, this candidate looks very good for HLA-A*11:01. 

BUT, there are two potentially problematic reference matches. First, a peptide sequence overlapping the predicted neoantigen peptide  is found in RASL10B. However in this case the register of the mutation in this matching sequence is shifted by two.  Presentation of that sequence is likely to be different (and was not found to be a strong binder).  Second, there is a longer peptide match to RHOT2 that contains the entirety of the neoantigen candidate peptide sequence (VVGARGVGK; including the R mutation itself).  RHOT2 is a RAS homolog family member.  According to GTEX, RHOT2 (and therefore the top peptide) is highly and ubiquitously expressed across many tissues. There are two other candidate strong binding peptides but they are just longer versions that contain much of the same sequence matching RHOT2 and therefore may have the same issue.

Perhaps the longer versions that do not match RHOT2 for the entire length are okay, but most of the sequence presented is still the same, and this sequence may have been presented in many tissues expressing RHOT2. 

My interpretation of this is that VVGARGVGK predicted to be expressed from mutant KRAS in the tumor cells and presented by HLA-A*11:01 matches exactly, by chance to VVGARGVGK from wild type RHOT2.  So presumably HLA-A*11:01 has been presenting that sequence from its expression in many non-tumor tissues and it would be subject to tolerance. The longer KRAS G12R peptides don’t match RHOT2 from end to end, but the portion facing the TCR seems like it would be much the same.  

GTEX expression data for RHOT:
![GTEX expression data for RHOT](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/f9990834-b85a-4cb3-a9ac-19e5b7bce242)

It is perhaps notable that in trial NCT03953235: “An improved V2 vaccine only targeted KRAS neoantigens (G12C/D/V, Q61H).”  Even the V1 vaccine design did not target G12R.

In a further twist of this candidate.  If we look at this position in RHOT2 there is a homozygous germline SNP that changes the homologous R to a C.  

![position in RHOT2 there is a homozygous germline SNP](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/fff57d29-bcb1-4724-acee-f8628363e86c)

So, this means even though there is a reference match when considering a wild type proteome that assumes the reference genome sequence, there is NOT a reference to this patient’s wildtype proteome.  So the candidate should be acceptable in the end. 

KRAS G12R peptide candidates and portion matching RHOT2:

CKVVGARGVGKSA RHOT2 sequence based on reference
  VVGARGVGK   (115 nM) 9-mer KRAS G12R candidate
 VVVGARGVGK   (164 nM) 10-mer KRAS G12R candidate
LVVVGARGVGK   (911 nM) 11-mer KRAS G12R candidate
CKVVGACGVGKSA RHOT2 sequence based on reference but with homozygous SNP applied



## pvacview screenshots -- different examples from different datasets
#### Leidos 5120-29 
SLC1A7 candidate (original status: Review) was dropped during the manual review. Short explanation: the Class I peptide which doesnt have reference match has bad binding affinity, Class II peptide has reference match. Long explanation: This candidate has 2 transcript sets (ENST00000620347.5 and ENST00000371494.9) , both transcript sets have a matched portion (LQALLIVL) with proteome reference. Class II peptide (QALLIVLATSSSSA, from ENST00000371494.9) has a complete overlap with the reference, thus was rejected. Class I peptide for HLA-C*03:04 (VLATSSSSATL, from ENST00000371494.9) can still be used (if we cut the 51 mer to exclude the reference match portion), however this peptide is a bad binder (median IC50 greater than 2,000 nM, with multiple algorithms reports high IC50). Class I peptide for HLA-A (ILQALLIVL from ENST00000620347.5) has good binding affinity but has a strong reference match.   

Further examples:
https://pvactools.readthedocs.io/en/latest/pvacview/pvacseq_module/pvacseq_vignette.html

