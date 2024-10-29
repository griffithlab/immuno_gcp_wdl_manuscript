# Troubleshooting README

The Immunotherapy Tumor Board (ITB) reviews all candidates and evaluates them as Accept, Reject, or Review. The candidates marked as Accept and Review will be reviewed in this 'Manual Review' section in two separate contexts, pVACview and IGV, to verify that they are good neoantigen candidates. 

During the manual review process, we will primarily use two files: the Annotated Neoantigen Candidates Excel sheet, the classI produced from the ITB review, and the Peptides 51mer sheet, which we will generate below. 
The Peptides 51mer sheet is generated in a certain format specified by the peptide manufacturer used for the neoantigen vaccine trials.

To generate the Peptides 51mer Excel sheet there are two steps: (1) generating the protein fasta and (2) generating the manual review files.

## TODO- NEED TO ADD SOME TROUBLESHOOTING EXAMPLES FOR FAILURES IN THE PIPELINE RUN ITSELF. (mentioning some general ideas here)

1. Most common errors related to user error when generating the YAML file. Check the cloud YAML file after cloudize workflows to make sure that all the file paths are in the google cloud bucket. And also check that all the files got uploaded. Also, check that the sample name in the readgroup is identical to the sample name in the field below.
2. Cromwell logging is difficult to parse, so there is not a single good way to troubleshoot a run failure. However, a good place to start is looking if there is a message pointing the user to a stderr log file from a step. 

## Preparing for ITB Review

After the pipeline run has been completed, we can verify that the outputs look normal and begin creating a summary of the pipeline to be presented at the ITB Review meeting where a panel of experts will review the predicted neoantigens. These meetings usually are about a 30 minutes to an hour per case so it is useful to create a case summary which we call a Genomics Review Report to orient everyone to the case. Some suggested information to include at the top of the report would be: an introduction to the case including things like cancer type, past treatments, any anynomallies and a list/links to the most important files used for this case (namely the IGV and pvacview files) .

### FDA Thresholds

Griffith lab discussions with the FDA have resulted in the creation of a set of values that set a very high standard of data quality. These metrics are pulled from several different files and then summarized in a table using the commands below:

```bash
mkdir $WORKING_BASE/../manual_review

docker run -it --env HOME --env WORKING_BASE -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts:latest /bin/bash

cd $WORKING_BASE/../manual_review

python3 /opt/scripts/get_FDA_thresholds.py -WB  $WORKING_BASE -f final_results
```

This is an example of how this data table should look for a very high data quality case. Note that the total reads often fail because the threshold is set at a very high level.

![FDA Quality Thresholds for Leidos-5120-28](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/0850f4d8-192f-4295-8d05-0a9907006ded)


This is an example of a case that suffered from low input tumor material resulting in tumor DNA coverage being deficient and uneven and suffering from a very high duplication rate. This likely contributed to a high false-positive rate with many variants not surviving basic filtering or manual review. Tumor RNA similarly suffered from an apparent high level of genomic contamination. As a result, several prioritized neoantigen candidates failed after dedicated IGV manual review, and the high duplication rate created situations where all variant support came from a set of identical/duplicate read alignments. 

![FDA Quality Thresholds for JLF-100-060](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/7fd55b36-586c-4882-b00e-cab71b8ea94a)

### Basic QC data review

Using similar steps from above we will also examine other basic QC data that we find most useful in evaluating the case.

```bash
mkdir $WORKING_BASE/../manual_review

docker run -it --env HOME --env WORKING_BASE -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts:latest /bin/bash

cd $WORKING_BASE/../manual_review

python3 /opt/scripts/get_neoantigen_qc.py -WB $WORKING_BASE -f final_results --yaml $WORKING_BASE/yamls/$CLOUD_YAML
```

We have created more concrete thresholds of evaluation for these metrics. 

Here is an example (case 5120-18) of a typical qc summary:

```
normal_dna_aligned_metrics.txt Unique Map Reads: 136,917,036 ( excellent )
tumor_dna_aligned_metrics.txt Unique Map Reads: 347,318,833 ( excellent )
tumor_rna_aligned_metrics.txt Unique Map Reads: 183,279,500 ( excellent )
normal_dna_aligned_metrics.txt Mapped Read Duplication Rate: 12.7410677726806 (%) ( very poor )
tumor_dna_aligned_metrics.txt Mapped Read Duplication Rate: 13.2917990956679 (%) ( very poor )
Relatedness: 0.994 ( excellent )
normal.VerifyBamId.selfSM Contamination: 0.00138 ( good )
tumor.VerifyBamId.selfSM Contamination: 0.00281 ( good )
The proportion of RNA reads mapping to cDNA sequence is 0.970317 (coding ( 0.896512 ) + UTR ( 0.073805 ) ( excellent )
trimmed_read_1strandness_check.txt: Data is likely RF/fr-firststrand
YAML file: immuno.strand: first
Total Number of somatic variants called: 67
```
#### The strand settings
A typical thing that the qc review can reveal is that the wrong strand setting was used in the yaml file. It is important to know the correct strand setting to determine the expression levels of gene, especially in cases of overlapping genes. So with an unstranded protocol, we do not map RNA reads to DNA by strand information. When the strand protocol is unknown, setting the strand setting in the yaml to unstranded is safest because it essentially does not consider the directionality of the reads during mapping. However, this does result in a loss of information. If you set your yaml to first strand when the protocol is second or unstranded, you could potentially lose whole genes depending on overlap...

**how to manually infer strandedness**

#### End Bias

Visually inspect the plot located at `qc/tumor_rna/rna_metrics.pdf`. If the RNA-seq data is of good quality from intact RNA the plot should have a “horseback” shape, representing lower coverage at the beginning and end of transcripts but quickly rising and remaining relatively high over the majority of the transcript positions. RNA-seq data made from highly degraded RNA combined with polyA selection or oligo-dT cDNA priming can have a heavily biased distribution instead. Such data can still produce gene expression estimates but will be unable to effectively verify the expression of somatic variant alleles.

An example of a good end bias plot:
![5120-30 end-bias horseshoe GOOD](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/2cea9675-5f30-45a2-b6ad-8afe458c7099)

An example of a poor end-bias plot. This case was run with the yaml set to "unstranded" when the detected strandedness of the RNA data was unstranded. **is this the reson of this plot???**
![jlf-100-067 end-bias BAD](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/649f9477-2d85-450c-8257-282c3be60abf)


### Fusion Review
Open the fusion inspector html '/gcp_immuno_workflow/rnaseq/fusioninspector_evidence/finspector.fusion_inspector_web.html', this web page will show possible fusions with evidence. A believable fusion would be one with more than 10 reads. Then we check to see if there are different left/right genes and different left-right chromosomes. Most of them are going to genes that are next to each other and it isread through rather than a real fusion. 

![fusion_table JLF-100-043](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/33f699ac-a555-4569-ac01-9797421a3f53)

The very first row is the only possible candidate with a lot of junctions reads but it's not spanning different chromosomes. This is an example of a read through:
![fusion_read_through JLF-100-043](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/cccd3c1a-41dd-4c28-9d40-0c9389fd72aa)

This example shows a possible frameshift in the first row, it has a ton of reads and fragments and is between two different chromosomes and genes.
![fusion table real jlf-100-056](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/c6e41e09-bd7a-4cf3-8228-3e7b3b45f0cd)

This is a real fusion, we see lots of spanning reads which are bringing the middle of the two genes together. 
![fusion_real jlf-100-056](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/3ddf57a4-2342-4750-91b4-9a8872e212d7)

**what do you do next? -- test for binding?? Run pvacfuse??**

#### JLF-100-037_mcdb041-original
Results from pVACfuse were used to include a peptide for a GFPT1::ENOX2 fusion in this tumor. The peptide sequence corresponds to the tumor-specific frameshift sequence created by the fusion event.

The ALK portion of the EML4::ALK fusion was extracted and used with pVACbind to nominate candidates for that driver event. In this case, while the ALK sequence is the portion presumably amplified/activated by the fusion event, the actual ALK sequence is NOT tumor specific.  It is simply the wild-type ALK sequence, from exon 20 to the end of the protein that is fused as the 3’ component of the EML4::ALK fusion. 

We also examined a fusion prediction for ZNF92::TDRD9 which had 44 junctions reads of support and looks real by manual review in fusion inspector and IGV. However this fusion is predicted to join the first exon of ZNF92 onto the second exon of TDRD9 and this is predicted to lead to an almost immediate stop codon.  The overall ORF would only be 10 AA or so and this is unlikely to be translated.

We also examined a fusion prediction for KDM5C::KMT2C which had 34 junction reads and 2 spanning reads of support.  It is predicted to lead to an almost immediate stop codon.  The neoORF segment would be: MLFHGCLS.  This truncation would be a drastic shortening of the normal KMT2C protein.  The only thing that is predicted to be a good binder is: MSHGVPMLF.  In other words only incorporating the first 3 AA of the novel sequencing arising from the fusion.  This is probably not worth targeting.

### HLA Allele Review
Make sure the HLA alleles that are being used in pVACview for the final selection of candidates match up with those expected for the case based on Optitype/PHLAT predictions from the data and clinical HLA typing results if available. Check whether the HLA alleles predicted for normal and tumor samples are in agreement.  Note any discrepancies.

Check the files in the location `gcp_immuno/final_results/hla_typing`. Make sure optitype_tumor_result.tsv and opitype_normal_result.tsv are the same, phlat_normal_HLA.sum and phlat_tumor_HLA.sum are the same. And all these calls match hla_calls.txt.

### Open case in pVACview

Take a look at the case in pVACview to make notes of anything that seems abnormal. Note how many candidates there are to review, the driving variants, framshift mutations, overall VAF/allele expression.

## ITB Review

Review the Genomic Review Report and note any concerning qc metrics, interesting candidates, or anything else abnormal about the case. Open pVACview and begin discussing candidates, generally taking note of:
- is it a good binder
- does it look well expressed
- is it a good class II binder?
- what are any potential probelms? TSL, reference match, problematic bases, etc
- Is in an anchor? Is it a stronger binder than wildtype?
- How many algorithms predict it a good binder?
- Are the elution scores good? (generally less that 1 percentile)
- Are there mutiple strong binding peptides arising from the variant?

Add comments and mark the evalution accordingly. Make sure to export your comments and evaluations after the ITB review is over. 

## Genomic Manual Review

### Genomic Manual Review Setup

#### Generate Protein Fasta

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

#### Generate the Peptide order form
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

### pVACview Review

Binding algorithm support, anchor position
Variants compared to CLE pipeline
IGV review
Somatic variant review (SOP)
Proximal variants review (also look for germline variants, neoantigen candidate may have to be changed to account for it)
RNAseq allele review
Transcript isoform review

##### Leidos 5120-29 
SLC1A7 candidate (original status: Review) was dropped during the manual review. Short explanation: the Class I peptide which doesnt have reference match has bad binding affinity, Class II peptide has reference match. Long explanation: This candidate has 2 transcript sets (ENST00000620347.5 and ENST00000371494.9) , both transcript sets have a matched portion (LQALLIVL) with proteome reference. Class II peptide (QALLIVLATSSSSA, from ENST00000371494.9) has a complete overlap with the reference, thus was rejected. Class I peptide for HLA-C*03:04 (VLATSSSSATL, from ENST00000371494.9) can still be used (if we cut the 51 mer to exclude the reference match portion), however this peptide is a bad binder (median IC50 greater than 2,000 nM, with multiple algorithms reports high IC50). Class I peptide for HLA-A (ILQALLIVL from ENST00000620347.5) has good binding affinity but has a strong reference match.  

#####  Assessing Read quality

the high duplication rate appeared to create situations where all variant support came from a set of identical/duplicate read alignments
![JLF-100-060 SLF2 variant high duplication rate](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/e5baa217-3cb0-4b4f-9ea4-f21dedee0133)

#####  Further examples:
https://pvactools.readthedocs.io/en/latest/pvacview/pvacseq_module/pvacseq_vignette.html

#####  General comment structure

When beginning review I used this general comment structure to make sure I was paying attention to the correct aspects.

1. Binding affinity for classI and %ile look good/ ClassII binding affinity and %ile
2. DNA, RNA VAFs and allele expression
3. Anchor position
4. Algorithm and elution 
5. Note anything else strange
   
### IGV Review

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

### Multiple Transcripts expressed

#### JLF-100-040_mcdb046
The CDKN2A gene has two supporting transcripts with slightly different 51-mer sequences so we have two transcripts in the final protein fasta list. Note that the 51-mer sheet does not include a classII sequence nor a candidate classII HLA allele for the transcript ENST00000579755.2 (those information are available for the best classII transcript which happens to be the other selected one).
	
BAZ2B also has two supporting transcripts. 

#### JLF-100-039_mcdb045 
The FAM76A has an unusual transcript that expresses the peptide in the intronic region and the used transcripts are not that one. 

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

### Alterations to best peptide that needs reanalysis of binding

#### JLF-100-039_mcdb045
The SLC9A6 on the other hand has a downstream missense variant (almost 4 bases away) that alters the expected best peptide, thus we do not have information of binding affinity, etc. for what would be the best new peptide. This candidate could potentially be rescued with additional effort if needed.

 #####  General comment structure

When beginning review I used this general comment structure to make sure I was paying attention to the correct aspects.

1. Somatic variant looks real (support in both RNA and DNA)
2. DNA, RNA VAFs and allele expression
3. Anchor position
4. Algorithm and elution 
5. Note anything else strange



