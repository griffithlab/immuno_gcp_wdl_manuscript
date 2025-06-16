# Manual Review REEADME

## Add something about running pVACvector

## Initial review of Immuno Pipeline Outputs
At this stage we mostly just want to confirm that the pipeline succeeded and that transferring the data files into the case folder on storage1 was successful.  Check the total number and total size of files.  In order to proceed these files must all be present: tumor/normal DNA BAMs, RNA BAM, variants VCF annotated for pVACseq, proximal variants VCF, annotated variants TSV, HLA typing results, FDA QC report files, pVACseq/pVACview results file. 
Presentation to the Immunogenomics Tumor Board (ITB)
Presentation of candidate neoantigens for selection by the Immunotherapy Tumor Board is performed using pVACview.  Most files needed to conduct this review are produced by the immuno pipeline and saved as results.  These consist of:
pVACview .R application files matched to the version of the results files
$sample-name.all_epitopes.aggregated.tsv (Class I)
$sample-name.all_epitopes.aggregated.metrics.json (Class I)
$sample-name.all_epitopes.aggregated.tsv (Class II)
Cancer Gene Census List in TSV format downloaded from Cosmic

Before the ITB meeting. Load these files in pVACview and make sure everything is working.  

During the meeting the “Evaluation” column will be used to mark each Pending candidate as Accept, Review, or Reject until sufficient candidates have been identified (some lower priority candidates may remain Pending). Once the case review is complete use the Export functionality to save the candidates and selections in both TSV and Excel format. These will be the ITB Reviewed Candidates that will be used during the following genomics review. 

ITB Evaluation Criteria
…
The following attempt to describe the criteria considered during immunogenomics tumor board meetings and enumerate any exceptions. 

Set Clonal DNA VAF. Based on the tumor DNA VAF of known drivers and/or consideration of the distribution of tumor variants you may choose to adjust the Clonal DNA VAF and “Recalculate Tiering” before beginning the review.  For example, if the tumor has a TP53 or KRAS driver variant with VAF of 25% and you believe it represents a heterozygous (non-CNV altered region) you might set this cutoff to that variant’s VAF.

Variant Type. Is the variant a SNV or an Indel. 

TSL (transcript support level). If the TSL is >1 or NA. Even if the candidate would otherwise be marked Accept it should be marked as Review and support for the specific transcript used to annotate predicted amino acid change should be evaluated using the tumor RNAseq data.

Pos (Position). Evaluate whether the variant is at an anchor site. If it is potentially at an anchor site, review the anchor heatmap view below.  If the variant is at an anchor, require that the wild type be a weak binder (>500 nm, ideally weaker than that). If the peptide is marked as Anchor Residue Fail in the “Transcript Set Detailed Data” view then mark the candidate as Reject.

Prob Pos (Problematic Positions). For peptide based studies, if any position in the core binding peptide has a Cysteine, then mark the candidate as Reject.

Num Passing Peptides. This is not an Accept/Reject criterion but a larger number of passing peptides is desirable, often the case for frameshift neoantigen candidates.

IC50 and %tile MT (median mutant peptide binding affinity prediction).  Should be less than 500nm or 1% for an Accept.  Exceptions are sometimes made: (a) if the variant is a known driver, (b) elution algorithms have a score < 1 %tile, (c) if there is disagreement between prediction algorithms and some predict that it is a strong binder (particularly NETMHCpan and MHCflurry).

RNA Expr (gene expression estimate). Should be >1 TPM for an Accept. A rare exception that might be evaluated during genomics review is when the RNA VAF and RNA coverage of the specific variant allele is strong but the gene expression was < 1 TPM.

RNA VAF (RNA variant allele fraction). Required >0 for an Accept.  Higher is better.  A low VAF may indicate sub-clonality or allele specific expression.  If the RNA depth is high, the DNA VAF is acceptable but the RNA VAF is low (e.g. <5%) then Reject (evidence of allele specific expression for the wild type allele).

Allele Expr (Allele Expression).  This metric is the product of the RNA VAF and Gene Expression value.  We have generally required a minimum Allele Expression value of 3 for Accept. Exceptions are sometimes made and a lower value is accepted in the 1-3 range, if the DNA VAF indicates that the variant is subclonal but otherwise a strong candidate. 

RNA Depth. This is not an Accept/Reject criterion but the value is used to interpret the RNA VAF.  For example, if RNA Depth < 10, the RNA VAF estimate is not very robust. Low RNA depth reduces confidence that the variant was actually detected. <DESCRIBE REVIEW SCENARIOS RELATING TO POOR RNA DEPTH>

DNA VAF. 

Tier.  

Transcript Sets. 

Transcript Biotype. For the most part, neoantigen predictions should correspond to transcripts that have been assigned the Biotype “protein_coding” by Ensembl. Additional biotypes may be considered (cautiously) if supported by review of the transcript annotation, RNA-seq data or other data (e.g. Mass Spectrometry) suggests the transcript actually does produce protein.  

Reference Matches.
Class II Binding Predictions.


Binding Affinity Algorithm Agreement. 


Elution Algorithm Scores.   


Cancer Gene Status.


Frameshift mutations.


Dinucleotide mutations.


Anchor Heatmap examination. 

Genomics Review (post ITB meeting)
Once the ITB meeting is complete and the preliminary set of candidates has been nominated, a genomics review of all candidates must be completed following the procedures and criteria laid out below.
Example Location of Files Needed for Genomics Review
Using case 5120-19 as an example:

/storage1/fs1/gillandersw/Active/Project_0001_Clinical_Trials/pancreas_leidos/analysis/TWJF-5120-19/gcp_immuno/final_results/ 

Variants Review Files
variants.final.annotated.tsv -> used to verify that Immuno variants were also called by CLE pipeline

pVACseq Review Files
pVACseq/mhc_i/5120_19_Tumor_FF.all_epitopes.aggregated.tsv -> Class I data
pVACseq/mhc_i/5120_19_Tumor_FF.all_epitopes.aggregated.metrics.json -> Class I data
pVACseq/mhc_ii/5120_19_Tumor_FF.all_epitopes.aggregated.tsv -> Class II data

IGV Review Files
annotated.expression.vcf.gz -> Immuno variant positions to load into IGV 
normal.cram -> Normal exome DNA alignments to load into IGV
tumor.cram -> Tumor exome DNA alignments to load into IGV
rnaseq/alignments/MarkedSorted.bam -> Tumor RNA-seq alignments to load into IGV
Genomics Review Checklist/Principles 
The following are high level descriptions of the kinds of detailed review that should be performed following selection of candidates by the ITB

Basic data QC review.  
Review relatedness, contamination, FDA QC reports, and purity metrics.  Note any red flags.

Summarize total number of uniquely mapping reads generated for Tumor/Normal exome and Tumor RNA-seq (i.e. not counting duplicates that arise from PCR amplification of the same source DNA fragment).
“Unique Mapped Reads” from file: “qc/fda_metrics/aligned_$sample_dna/table_metrics/$sample_dna_aligned_metrics.txt”
Qualitative description of these counts: 
<40,000,000 is poor
40,000,000 - 50,000,000 is acceptable
50,000,000 - 100,000,000 is good
> 100,000,000 is excellent
Summarize duplication rates for tumor/normal DNA samples
“Mapped Read Duplication” from file: “qc/fda_metrics/aligned_normal_dna/table_metrics/normal_dna_aligned_metrics.txt” 
Qualitative description of these rates:
>75% is very poor
50-75% is poor
30-50% is acceptable
20-30% good
<20% is excellent

Check Somalier results for sample tumor/normal sample relatedness.
“Relatedness” column from file: “qc/concordance/concordance.somalier.pairs.tsv” 
Qualitative description of these rates:
>97.5% is excellent
95%-97.5% is good
90-95% is concerning
<90% is very concerning

Check VerifyBamID results for contamination of both tumor and normal samples
“FREEMIX” column from file (column 7): “qc/$sample/$sample.VerifyBamId.selfSM”
Qualitative description of these rates:
<0.025 is good
0.025-0.05 is concerning
>0.05 is very concerning

Check RNA-seq metrics for % reads aligning to transcripts
Find the sum of “PCT_CODING_BASES” and “PCT_UTR_BASES” from file: qc/tumor_rna/rna_metrics.txt
Helpful command: “cut -f 17,18 rna_metrics.txt” 
>90% excellent
75-90% good
50-75% acceptable
<50% sub-optimal

Check for evidence of end bias in the aligned RNA seq data using file: “qc/tumor_rna/rna_metrics.pdf”
Visually inspect this plot. If the RNA-seq data is of good quality from intact RNA the plot should have a “horse back” shape, representing lower coverage at the beginning and end of transcripts but quickly rising and remaining relatively high over the majority of the transcript positions. RNA-seq data made from highly degraded RNA combined with polyA selection or oligo-dT cDNA priming can have a heavily biased distribution instead. Such data can still produce gene expression estimates but will be unable to effectively verify expression of somatic variant alleles.

Check that the correct RNA strand setting was used in the pipeline YAML file. For example, the 
Detected strand file: qc/tumor_rna/trimmed_read_1strandness_check.txt
YAML file: workflow_artifacts/$case_immuno_cloud-WDL.yaml (grep for “strand”)

Summarize total variants called and total neoantigen variants called (i.e. the subset of variants that could lead to neoantigens).  Briefly, total variants is how many somatic mutations were detected in the tumor. Total neoantigen variants is the subset of these that lead to possible neoantigen candidates (i.e. those that cause protein coding changes in known genes). Note that the number of neoantigen candidates selected by the Immunotherapy Tumor Board will be a small subset of this number.
Total variants called can be obtained from “variants.final.annotated.tsv” (previously “pvacseq.annotated.tsv”). Simply the number of rows in this file minus the header.
Total neoantigen variants called can be obtained from pVACview (total candidates) or by looking at the pVACseq aggregated report file to get this number.

FDA Quality Thresholds
The following subset of QC values can be obtained from the FDA report tables generated by the pipeline and used to determine whether the case meets basic data quality criteria as described in documentation provided to the FDA.

Criteria
Threshold for pass/failure
TOTAL_READS (tumor DNA)
>250M
TOTAL_READS (normal DNA)
>100M
TOTAL_READS (tumor RNA)
>200M
PCT_PF_READS_ALIGNED (tumor/normal DNA)
>95%
PCT_PF_READS_ALIGNED (tumor RNA)
>85%
PCT_USABLE_BASES_ON_TARGET  (tumor/normal DNA)
>20%
PCT_EXC_OFF_TARGET  (tumor/normal DNA)
<60%
PERCENT_DUPLICATION  (tumor/normal DNA)
<40%
MEAN_TARGET_COVERAGE (tumor DNA)
>250x
MEAN_TARGET_COVERAGE (normal DNA)
>100x
PCT_TARGET_BASES_20X  (tumor/normal DNA)
>95%
PCT_READS_ALIGNED_IN_PAIRS  (tumor/normal DNA)
>95%
MEAN_INSERT_SIZE (tumor/normal DNA)
125 - 300bp
PF_MISMATCH_RATE_1 (tumor/normal DNA)
<0.75%
PF_MISMATCH_RATE_2 (tumor/normal DNA)
<1.00%
Genotype Concordance (DNA pair)
>95%
Contamination Estimate (tumor/normal DNA)
<7.5%



Tumor type / driver variant review
Given the reported tumor type do we see evidence for expected somatic driver variants?  Enumerate these in the report and if there are no plausible driver variants, raise this as a possible red flag.  This review can be conducted in pVACview with CancerGeneCensus genes loaded.

Other strategies for assessing drivers:
Use the variants TSV file and select all variant types that would directly impact protein sequence
Look for variants previous observed in Cosmic or pathogenic according to ClinVar
Intersect the genes with these mutations with Cancer Gene Census
Take the candidate genes that result from the previous two queries and search for these in CBioPortal after selecting studies that match the cancer type of the patient

HLA allele review for sample/data mixup
Make sure the HLA alleles as being used in pVACview for the final selection of candidate match up with those expected for the case based on Optitype/PHLAT predictions from the data and clinical HLA typing results if available.

HLA allele review of normal vs tumor
Check whether the HLA alleles predicted for normal and tumor samples are in agreement.  Note any discrepancies.

Review multi-algorithm support for strong binding affinity
We are often using the “lowest” score method to prioritize candidates (the other, more conservative option would be “median”).  This means the peptide with the lowest IC50 from *any* algorithm. This leads to more candidates but can lead to the situation where a peptide is nominated and most algorithms actually predict that it is a poor binder and a minority of algorithms (perhaps only one) predict that it is a strong binder.  In those cases, is the candidate still acceptable? Is there another peptide where the best score was slightly higher but there is better agreement across algorithms?
MHCnuggets has often been observed to be an outlier.  If only this algorithm suggests a peptide is a strong binder, these candidates should probably be failed in review.
NetMHCPan and MHCFlurry have been reported by some benchmarking exercises as among the better performing algorithms
In cases where binding algorithms have high disagreement, additional weight should be given to presentation algorithms trained on peptide-elution mass spectrometry data

Compare variant results to CLE pipeline
For all variants with neoantigen candidates, check if the variant was also called by the CLE pipeline. This status is automatically provided in the “VALIDATED” column in variants.final.annotated.tsv file produced by the immuno pipeline. 

IGV Review
You SMB mount to obtain access to the case data files on your local computer and create an IGV session with the following 6 components in order: 
Final somatic variant VCF from immuno pipeline (annotated.expression.vcf.gz)
CLE variant VCF (annotated_filtered.vcf.gz.commented.gz.commented.gz)
Normal exome DNA alignments (normal.cram)
Tumor exome DNA alignments (tumor.cram)
Tumor RNA-seq alignments (rnaseq/alignments/MarkedSorted.bam)
Ensembl v105 transcript annotations (Homo_sapiens.GRCh38.105.sorted.coding.gtf) 

Save the session file on storage1 for others to use.  Use this session to perform the following specific review activities and make note of the findings in the candidate review spreadsheet.

Somatic variant review. For any peptide candidate being considered for inclusion in the design, manually review the underlying somatic variant in IGV to make sure it appears to be a real somatic variant.  Our published guidelines can be used for this step. 
Proximal variants review. Check for any nearby variants that might impact prediction of the peptide sequence. Most of these should be handled automatically by pVACview if it was run with the proximal variants option. The definition of “nearby” for the purposes of proximal variants review is any variant capable of impacting either the predicted class I epitope sequence (8-11 amino acids), the class II epitope sequence (12-18 amino acids) or the long peptide sequence that will be used for DNA/RNA vector or peptide sequencing (35-50 amino acids). If a variant is NOT near the edge of an exon, the proximal space to review may be up to 75 bp in either direction (potentially longer for a frameshift variant). If the variant is near the edge of an exon, the proximal space to review can be much larger as adjacent exons may be separated by large introns.  
RNA-seq allele review. Check for expression of the allele.  If the candidate was selected, but had poor read support for the mutant allele, check whether this could be explained by a local coverage or alignment issue, where the gene is otherwise well expressed.
Transcript isoform review. Check if the RNA transcript the candidate peptide is extracted from is actually expressed. Does the RNA-seq data support the expected nearby exon-exon junctions?  Is there evidence for alternative splicing that might influence expression of the mutant allele or the structure of the transcript that expresses (possibly impacting the presumed peptide sequence)? When performing this review the following IGV display settings are useful: (1) “View as pairs” and (2) “Color alignments by” -> “first-of-pair strand”. Note that for this step it will be necessary to load the correct Ensembl GTF of transcripts in IGV. A sorted and indexed GTF (v105; the version currently used in the pipeline) has been shared in Google Drive, Google Cloud Bucket and on Storage1: /storage1/fs1/gillandersw/Active/Project_0001_Clinical_Trials/annotation_files_for_review/. 

Long Peptide Extraction and Annotation
Using the selected transcript, extract 51-mer peptide sequences that contain the candidate neoantigen.  Annotate this sequence with “best” classI and classII binding peptides.  The selection of long peptide sequences chosen should reflect the final conclusion of the ITB review and genomics review.
Use “pvacseq generate_protein_fasta” to create a fasta file that has the long peptide sequences needed.
For convenience you can create two versions of the Fasta file:
1. One that has a mutant peptide sequence for every mutation that was evaluated “Accept,Review” during the ITB meeting and based only on the top transcript. This should be close to the final set of sequences you want. 
2. One that has a mutant peptide sequence for every mutation, generated for every transcript.  This may be needed to manually adjust the final set of peptides to account for transcript expression, alternative splicing, etc.

Final Reporting
All of the above findings should be summarized in a report and shared with the ITB group for final review prior to ordering the vaccine.  The report files should be uploaded to the appropriate patient folder in Wustl Box and consist of the following five components:
Peptide Fasta. A plain text fasta file of each peptide sequence. This will be used for final blast checks and can be used with Pepstats to obtain molecular weights for each peptide for the ordering form. 
Note that if any changes to peptides are made manually during the genomics review (e.g. to select a different reference transcript, or to correct for a complex variant or proximal variant), this file must be updated to have a one-to-one relationship with the Long Peptides Spreadsheet.
The starting point for this file is the “annotated_filtered.vcf-pass-51mer.fa” file produced by the “pvacseq generate_protein_fasta” tool described above. Example location: $case/generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa

Long Peptides Spreadsheet. An annotated peptide sequence file used to create the final vaccine manufacturing order form. The starting point for this spreadsheet is the “.fa.manufacturability.tsv” produced by the “pvacseq generate_protein_fasta” tool described above. Example location: $case/generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv.     The Long Peptide Spreadsheet should provide the following:
Peptide highlights. ClassI (red text), ClassII (bold text) and mutant peptides (underlined) and if using a long peptide vaccine Cysteines (C) indicated with a larger font. Please note the underline convention for these variant types:
SNVs - underline the single mutant amino acid
In-frame insertions - underline the inserted amino acids
In-frame deletions - underline the one amino acid to the left and one amino to the right of the deleted amino acid position
Frameshift Insertions or Deletions - underline the entire novel amino acid sequence (i.e. from where the sequence starts to diverge from wild type until the end of the peptide sequence).
In-frame RNA fusions - underline one amino acid on the left side of the fusion junction and one amino acid on the right side of the fusion junction.
Frameshift RNA fusions - underline the entire novel amino acid sequence (i.e. from where the sequence starts to diverge from wild type until the end of the peptide sequence).
Top class I/II HLA alleles
List the best class I and class II allele separated by a “/”. List the class I allele if median affinity < 1000 nm OR percentile < 2%. List the class II allele if percentile < 2%. Note that the class I/II peptide should only be red/bold if it meets these criteria.
Molecular weight
Get the molecular weight for each long peptide sequence by using the EMBOSS Pepstats program.

Reviewed Candidates Table (MHC Class I). A spreadsheet with the reviewed class I candidates from pVACview along with any notes on each candidate. Each row corresponds to a single variant leading to potential neoantigens, with the top candidate provided.  Each should include the final evaluation: Reject or Accept.  Each should also include any detailed notes from the pVACview and IGV reviews that document the rationale for Accepting or Rejecting the candidate.
The starting point for this spreadsheet is the TSV that is exported from pVACview at the end of the review in the Immunotherapy Tumor Board meeting. Example storage1 location and file naming convention: $case/itb-review-files/mcdb047.revd.Annotated.Neoantigen_Candidates.xlsx

Reviewed Candidates Table (MHC Class II). A spreadsheet with the best class II candidate peptides. Each row corresponds to a single variant leading to potential neoantigens, with the top candidate provided. 
The starting point for this spreadsheet is the class II aggregate report from pVACseq that came from the WDL pipeline (or a manual run of pVACseq if that was required). Example storage1 location: mcdb046/gcp_immuno/final_results/pVACseq/mhc_ii/mcdb046-tumor-exome.all_epitopes.aggregated.tsv
 
Executive Report. A narrative report with an executive summary of key findings from each of the major genomics review steps above. This should be created as each step is completed. In addition to summarizing each step, note any special situations with individual candidates.  For example, if a complex variant was manually corrected, if candidates were altered or added to reflect transcript annotation issues or observed alternative splicing, special consideration given to driver variants, candidates selected based on class II alone, reference matches, etc.  If any peptides are determined or corrected by manual work or ad hoc analysis, describe those in detail here along with visualizations if helpful. 

Communication
Email the ITB team to let everyone know that the post-ITB genomics review is complete. 
