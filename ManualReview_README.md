# Immunogenomics Review 

The purpose of the immmunogenomics review is to confirm that the results of the computational pipline are correct and reasonable. This process involves checking if the output contains the correct files, reviewing quality control (qc) metrics, creating a written report, presenting results to a panel of experts, and using visualization tools, IGV and pVACview, to manually verify the neoantigen candidates.These steps assure that the neantigen candidates of interest are of the highest confidence.

## Initial review of Immuno Pipeline Outputs

At this stage we mostly just want to confirm that the pipeline succeeded and that transferring the
data files into the case folder on storage 1 was successful. Check the total number and total
size of files. In order to proceed, these files must all be present: tumor/normal DNA BAMs, RNA
BAM, variants VCF annotated for pVACseq, proximal variants VCF, annotated variants TSV,
HLA typing results, FDA QC report files, pVACseq/pVACview results file.

## Presentation to the Immunogenomics Tumor Board (ITB)

Presentation of candidate neoantigens for selection by the Immunotherapy Tumor Board is
performed using pVACview. Most files needed to conduct this review are produced by the
immuno pipeline and saved as results. These consist of:
- pVACview .R application files matched to the version of the results files
- $sample-name.all_epitopes.aggregated.tsv (Class I)
- $sample-name.all_epitopes.aggregated.metrics.json (Class I)
- $sample-name.all_epitopes.aggregated.tsv (Class II)
- Cancer Gene Census List in TSV format downloaded from Cosmic

Before the ITB meeting. Load these files in pVACview and make sure everything is working.

During the meeting the “Evaluation” column will be used to mark each Pending candidate as
Accept, Review, or Reject until sufficient candidates have been identified (some lower priority
candidates may remain Pending). Once the case review is complete use the Export functionality
to save the candidates and selections in both TSV and Excel format. These will be the **ITB
Reviewed Candidates** that will be used during the following genomics review.

## ITB Evaluation Criteria

The following sections describe the criteria considered during immunogenomics tumor board
meetings, provide examples, define Accept/Reject criteria where possible and enumerate
special cases and exceptions.

**Set Clonal DNA VAF**

Based on the tumor DNA VAF of known drivers and/or consideration of the distribution of tumor
variants you may choose to adjust the Clonal DNA VAF and “Recalculate Tiering” before
beginning the review. For example, if the tumor has a TP 53 or KRAS driver variant with VAF of
25 % and you believe it represents a heterozygous (non-CNV altered region) you might set this
cutoff to 25 %.

**Variant Type**

Consider whether the variant is a SNV, InDel or DNV. Variants that are insertions, deletions or
nucleotide variants have different interpretations from the perspective of DNA and RNA VAF,
reference matches, number of core binding peptides expected, etc.


[**TSL (transcript support level)**](https://grch37.ensembl.org/info/genome/genebuild/transcript_quality_tags.html)

Take note if the TSL is > 1 or NA. Even if the candidate would otherwise be marked Accept it
should be marked as **Review** and support for the specific transcript used to annotate the
predicted amino acid change should be evaluated using the tumor RNAseq data. TSL values
greater than 1 or NA are often acceptable and simply reflect incomplete annotation of the
transcripts of a gene (an ongoing, partially manual, curation process). It can be helpful to refer
to the “transcript table” in Ensembl where additional annotation tiering information is available.
For example, it is possible for a transcript to have a TSL > 1 , but for that transcript to also be
categorized as MANE Select, GENCODE Primary or Ensembl Canonical, any one of which
would increase confidence in the quality of the transcript annotation. A screenshot of an
example Ensembl transcript table is provided below for reference:

<img width="1399" alt="ensemble transcript" src="https://github.com/user-attachments/assets/b3535886-6c38-4f73-880a-d9b021648d1f" />


**Pos (Position)**

Evaluate whether the variant is at an anchor site. If it is potentially at an anchor site, review the
anchor heatmap view below. If the variant is at an anchor, require that the wild type be a weak
binder (> 500 nm, ideally weaker than that). If the peptide is marked as Anchor Residue Fail in
the “Transcript Set Detailed Data” view (the peptide and HLA listed under Anchor Residue Fail
column will be highlighted red) then mark the candidate as **Reject**.

**Prob Pos (Problematic Positions)**

In past testing it has been observed that peptides containing cysteines may result in difficulties
during synthesis, disulfide bond formation, or reduced stability during storage. Peptides
containing cysteines may later fail identity/purity tests (e.g. by LC MS) during drug product
testing (which may be requested by the FDA).

For most studies using a synthetic long peptide approach for vaccine delivery, if any position in
the core binding peptide has a cysteine, then mark the candidate as **Reject**. Some long peptide
vaccine studies allow some amount of cysteines in the synthesized peptide (e.g. 1 - 2 in the long
peptide sequence) and may attempt to mitigate the potential issues above (e.g by strategic
pooling, such as avoiding multiple cysteine containing peptides in the same pool).

For studies using DNA/RNA vectors for vaccine delivery, no problematic amino acids will be
defined in the pVACtools analysis and this criteria can be ignored.


**Num Passing Peptides**

This is not an **Accept** / **Reject** criterion but a larger number of passing peptides is desirable, as is
often the case for frameshift neoantigen candidates. Note that the number of peptides reported
as Passing in pVACview and aggregated reports depends on the thresholds selected. Relevant
pVACtools parameters include: _- -binding-threshold_ , _- -percentile-threshold_ ,

_- -allele-specific-binding-thresholds_ and _- -top-score-metric_.

**IC 50 and %tile MT (median mutant peptide binding affinity prediction)**

The median IC 50 should be less than 500 nm or the median percentile less than 1 % for an
Accept. Exceptions are sometimes made if: (a) the variant is a known driver, (b) elution
algorithms (e.g. BigMHC_EL, MHCflurryEL Presentation, NetMHCpanEL) have a score with <
1.0 %tile, (c) if there is disagreement between prediction algorithms and some predict that it is a
strong binder (particularly NETMHCpan and MHCflurry as these algorithms have performed well
in benchmarking exercises).

**RNA Expr (gene expression estimate)**

The gene expression estimate should be > 1 TPM for an Accept. A rare exception that might be
evaluated during genomics review is when the RNA VAF and RNA coverage of the specific
variant allele is strong but the gene expression was < 1 TPM.

**RNA VAF (RNA variant allele fraction)**

An RNA VAF > 0 is generally required for an Accept. Higher is better. A low VAF may indicate
sub-clonality or allele specific expression. If the RNA depth is high, the DNA VAF is acceptable
but the RNA VAF is low (e.g. < 5 %) then Reject (evidence of allele specific expression for the
wild type allele). In rare cases, exceptions may be made if: (a) the variant is a known driver, (b)
the gene is expressed at a high level; and (c) the low RNA VAF may be explained by lower
coverage for technical reasons (e.g., in a region suffering from end bias).

**Allele Expr (Allele Expression)**

This metric is the product of the RNA VAF and Gene Expression value. We have generally
required a minimum Allele Expression value of 3 for Accept. Exceptions are sometimes made
and a lower value is accepted in the 1 - 3 range, if the DNA VAF indicates that the variant is
subclonal but otherwise a strong candidate.

**RNA Depth**

This is not an Accept/Reject criterion but the value is used to interpret the RNA VAF. For
example, if RNA Depth < 10 , the RNA VAF estimate is not very robust. Low RNA depth reduces
confidence that the variant was actually detected.

If the total number of variant allele supporting RNA reads is small (e.g. < 5 ), IGV review should
be used to carefully assess the quality of this small number of reads. The number of unique
cDNA fragments supporting the variant should be examined (e.g. if the supporting RNA count is
2 , it is possible that there is just 1 unique cDNA fragment, sequenced in both directions).


In the case of small numbers of supporting RNA reads, the alignments of these reads should be
examined to confirm that they are consistent with a mature RNA sequence. For example, if the
RNA reads containing the variant are spliced across introns (exon-exon junctions) this would
represent stronger support than if they are contained within an exon. Similarly, if the RNA-seq
data was generated in a strand-specific way, RNA reads containing the variant that also match
the expected expression strand would represent stronger support than if they have the incorrect
strand. If RNA reads containing the variant are not correctly spliced (e.g. span into the intron
instead of connecting to the next exon) these should not be considered supporting the variant
even though they contributed to a non-zero VAF.

**DNA VAF**

By definition DNA VAF must be > 0 if a variant was called. The accuracy of a DNA VAF estimate
depends on the depth of reads covering the variant position.

A higher DNA VAF is preferable because it suggests higher tumor clonality. A higher DNA VAF
can also result if the variant allele has been amplified or is affected by loss-of-heterozygosity
(LOH). For example, if a heterozygous variant arises in a tumor and the chromosome (or a
region of it) that harbors the variant is duplicated, this can result in higher DNA VAF. This
scenario is also desirable from the perspective of neoantigen prioritization.

A lower DNA VAF is less desirable because it could indicate a sub-clonal variant (i.e. one that is
not present in all the cells of the tumor). Low DNA VAFs are also more likely to be false positive
somatic variants, although this is dependent on the absolute depth of sequencing achieved. If
the sequence error rate of the platform is low, and the depth is high, somatic variants can be
robustly identified down to low VAFs (e.g. 1 - 2 %). A default conservative threshold for allowing
somatic variant calls is often set at 5 % but may be lowered to 1. 5 - 2. 5 % in cases where tumor
purity is low.

Correct interpretation of individual DNA VAFs generally must consider multiple factors including:
(a) purity of the tumor (aka tumor cellularity); (b) ploidy of the tumor (e.g. copy number
gains/losses that may influence VAFs); (c) clonality of the variant; (d) random sampling error. A
detailed discussion of these interpretations is beyond the scope of this document but many tools
and publications are available.

**Tier**

Each neoantigen candidate variant is automatically placed into one of several tiers based on a
simple heuristic: Pass, LowExpr, NoExpr, SubClonal, Poor, Anchor. The Tier is evaluated for the
top peptide arising from each variant.

If the Tier is Pass, the variant can generally be accepted with limited review. However, there are
a few issues that may require even a Pass candidate to be rejected including: (a) problematic
amino acids (e.g. Cysteines, if defined by the user) in the core binding peptide; (b) problematic
amino acids that are outside the core binding sequence but closely flanking it such that a long
peptide sequence is not possible (c) reference matches, (d) other issues that come up during genomic review, 
including that the somatic variant giving rise to the neoantigen is a false
positive.

The additional tiers (LowExpr, NoExpr, SubClonal, Poor, Anchor) are designed to draw attention
to specific areas of consideration during the ITB review or genomics review.

The Poor tier indicates that either the candidate does not meet basic binding/presentation
criteria or it does but has multiple other issues.

**Transcript Sets**

The recommended approach for pVACtools analysis is to annotate somatic variants with all
possible transcript annotations according to VEP (using the - -flag-pick option to note the top
transcript consequence according to VEP but still retain the others). The result of this approach
is that neoantigen predictions for a single variant can vary from one transcript annotation to
another. A neoantigen may arise in an exon of one transcript (where the variant is considered
missense) but be absent in another transcript if the exon is skipped (where the variant is
considered intronic). Furthermore, a variant with a consistent predicted amino acid change in
two transcripts may still lead to distinct short (e.g. 8 - 11 - mer) or long (e.g. 25 - 35 - mer) peptides if
the variant is close to the edge of an exon and the RNA splicing of downstream/upstream exons
varies. These are just two examples of the many impacts that alternative transcripts may have
on neoantigen identification and prioritization. While not common, it is possible for entirely
distinct neoantigens to arise from a single somatic variant when multiple RNA transcripts are
simultaneously expressed.

In pVACview, a “transcript set” defines the collection of transcripts that give rise to an identical
list of predicted peptides. If more than one transcript set exists, the set producing the top single
peptide candidate is automatically selected. Two transcript sets can have overlapping peptides
but the overall complement of peptides from each set must differ. I.e. there must be at least one
peptide sequence that is only possible in the context of one transcript set and not the other.

When multiple transcript sets are observed for a somatic variant, they should each be
considered to determine whether multiple long peptides could be targeted in the neoantigen
therapy. Alternative candidate peptides arising from each set must meet the usual criteria for
binding, presentation, expression, transcript structure support, etc.

As a final note, all of the above functionality and guidance on alternative transcripts applies only
to known alternative RNA transcripts. Somatic variants that alter RNA splicing in a tumor
specific fashion can also be a rich source of neoantigen but must be analyzed by a dedicated
approach (e.g. using pVACsplice).

**Transcript Biotype**

For the most part, neoantigen predictions should correspond to transcripts that have been
assigned the Biotype “protein_coding” by Ensembl. Additional biotypes (e.g. nonsense mediated
decay) may result in neoantigen peptide predictions and may be considered (rarely and


cautiously) if review of the transcript annotation, RNA-seq data or other data (e.g. Mass
Spectrometry) suggests the transcript actually does produce protein. It is often helpful to review
the “transcript table” in Ensembl (see example above in the TSL section).

**Reference Matches**

Reference Matches are an Accept/Reject criterion. Reference matches refer to one or more
sequence matches between the candidate neoantigen peptide and other regions of the
reference proteome (e.g. we typically compare to all known Ensembl protein sequences). The
reference match analysis performed by pVACtools first attempts to identify the expected match
of the mutant sequence to its corresponding wild type gene location. This result is excluded
from the Reference Match result and any remaining matches to other genes or other positions
within the source gene are reported.

Generally a candidate with reference matches is marked as Reject. A larger number of
reference matches is considered more problematic. A large number of reference matches (e.g.
\> 3 ) is often observed when the neoantigen peptide has low complexity (e.g. contains a
homopolymer stretch) or when it arises from a large gene family with a high degree of homology
across family members.

Correct interpretation of Reference Matches in pVACview requires understanding of three
related amino acid sequences and how they relate to each other: the Best Peptide, the Query
and the Matched Peptide.

The Best Peptide sequence is the sequence that was evaluated as a neoantigen candidate (e.g.
an 8 - 11 - mer mutant peptide for class I neoantigens. The Query Sequence is used to search for
reference matches related to this sequence. This sequence is created by identifying the position
of the mutation and adding a flanking sequence of 7 amino acids evenly to either side. The
Matched Peptide is a substring of the Query Sequence. The Matched Peptide comes from an
unexpected location in the reference proteome (i.e anywhere other than the region the
neoantigen Best Peptide sequence comes from). The Matched Peptide may have complete or
partial overlap with the Best Peptide sequence but the overlap must contain the mutant amino
acid position(s). Refer to the examples below to visualize the relationship between the Best
Peptide, the Query and the Matched Peptide.

In certain circumstances a candidate with reference matches may be considered acceptable.
Generally a candidate with a reference match will only be acceptable if the number of reference
matches is low (e.g. 1 ) and the degree of overlap between the Best Peptide candidate and the
Matched Peptide is limited (e.g. less than 50 % of amino acids of the Best Peptide occur in the
Matched Peptide).

In the following examples, the Best Peptide sequence is in bold, the mutant amino acid is
marked in red and the substring that matches between all three sequences is underlined.

**Example 1** (from hcc 1395 benchmark analysis):
<img width="951" alt="Ref Match SORBS3" src="https://github.com/user-attachments/assets/148197fd-2ce2-465d-a263-cdb951d486c4" />


**Example 2** (from hcc 1395 benchmark analysis):
<img width="971" alt="REf Match FAM217B" src="https://github.com/user-attachments/assets/3f63b2b1-b787-4484-aa52-39ea14ec38d5" />

The general principle of Reference Match analysis is that we want to avoid the situation where
the neoepitope sequence presented to the T cell by the Best Peptide could also be presented
by the reference peptide. Interpretation of this concept is often aided by referring to the Anchor
Heatmap view in pVACview (see below for more details).

A known complication of the reference match analysis is that for inframe insertions and
deletions, the mutant peptide will sometimes be reported as having a reference match to the
same gene it derives from. This result occurs because of the ambiguity in interpreting the
location of insertions and deletions, particularly where there is similarity between the
inserted/deleted amino acids and the surrounding reference amino acids.

Refer to the pVACtools and pVACview documentation and tutorials for more information on the
determination and interpretation of reference matches.

**Class II Binding Predictions**

In the case where a neoantigen candidate does not have any qualifying class I peptide
predictions (no convincing MHC binding, presentation and immunogenicity), a candidate may
still be Accepted based on class II MHC predictions. Relaxed criteria to Accept a candidate
based on class II binding are: < 500 nm OR percentile < 2 %. Class II predictions tend to result
in systematically lower binding affinity values for many class II alleles and a more conservative
criteria requires the percentile to be < 2 %.


*Note on class II pairing notations*

In class I MHC notation, a complex is formed between each class I HLA molecule and B 2 M (e.g.
HLA-B* 45 : 01 - B 2 M). However, since the involvement of B 2 M is a constant its presence is
implied and it is not listed in the HLA notation. In class I, the MHC peptide binding groove is
internal to the protein encoded by the HLA gene.

In contrast to class I, all class II MHC molecules operate as a dimer with the MHC peptide
binding groove formed by the dimerization. Each class II complex is described as an alpha-beta
pair. The three HLA genes of primary interest for class II binding and presentation are HLA-DP,
HLA-DQ, and HLA-DR.

**HLA-DP** and **HLA-DQ** complexes have a dimer notation, explicitly listing the alpha-beta pair of
each complex. Valid HLA-DP and HLA-DQ dimer notations for neoantigen analysis therefore
take the form:

[HLA-DPA 1 *NN:NN - HLA-DPB 1 *NN:NN]

[HLA-DQA 1 *NN:NN - HLA-DQB 1 *NN:NN].

An individual typically has up to two alleles for DPA 1 , DPB 1 , DQA 1 , and DQB 1. All DPA 1
alleles can pair with all DPB 1 alleles. All DQA 1 alleles can pair with all DQB 1 alleles.

Note that we have encountered algorithms that will report binding predictions for single
components of HLA-DP and HLA-DQ genes (only alpha or only beta). Such predictions are
difficult to interpret and should be ignored.

Example HLA types for the HLA-DP gene:
Patient 1. DPA 1 * 01 : 03 , DPA 1 * 02 : 01 DPB 1 * 04 : 02 , DPB 1 * 14 : 01

Valid pairings:
DPA 1 * 01 : 03 - DPB 1 * 04 : 02

DPA 1 * 01 : 03 - DPB 1 * 14 : 01

DPA 1 * 02 : 01 - DPB 1 * 04 : 02
DPA 1 * 02 : 01 - DPB 1 * 14 : 01

Example HLA types for the HLA-DQ gene:
Patient 1. DQA 1 * 01 : 01 , DQA 1 * 01 : 02 DQB 1 * 05 : 01 , DQB 1 * 06 : 02
Valid pairings:
DQA 1 * 01 : 01 - DQB 1 * 05 : 01

DQA 1 * 01 : 01 - DQB 1 * 06 : 02

DQA 1 * 01 : 02 - DQB 1 * 05 : 01
DQA 1 * 01 : 02 - DQB 1 * 06 : 02


In the case of **HLA-DR** genes, by convention, only the beta component is listed because
HLA-DR alpha is not functionally variable across individuals (i.e. the HLA-DRA component is
constant/implicit). Valid HLA-DR dimer notations for neoantigen analysis therefore take the form:
DRB 1 *NN:NN, DRB 3 *NN:NN, DRB 4 *NN:NN and DRB 5 *NN:NN. An individual typically has up
to two alleles for the DRB 1 gene and up to two alleles drawn collectively from the DRB 3 / 4 / 5
genes.

Example HLA types for the DR genes:
Patient 1. DRB 1 * 07 : 01 , DRB 1 * 14 : 07 DRB 3 * 02 : 02 , DRB 4 * 01 : 03

Patient 2. DRB 1 * 03 : 01 , DRB 1 * 15 : 01 DRB 3 * 01 : 01 , DRB 5 * 01 : 01
Patient 3. DRB 1 * 01 : 01 , DRB 1 * 15 : 01 DRB 5 * 01 : 01

**Types of peptide:MHC prediction scores**

Several types of prediction scores have been described as relevant to pMHC prioritization.
Included among these are distinct concepts: binding, stability, cleavage, transport, processing
and immunogenicity. A description of the rationale behind each of these types of algorithms is
far beyond the scope of this document, but in very simple terms:
- **Binding** (aka binding affinity) algorithms (e.g. NetMHCpan, MHCflurry) directly model the
binding of a short peptide sequence to the class I/II MHC groove. These algorithms are
typically trained on in vitro binding data generated by testing individual peptides (from
peptide libraries) against individual MHC molecules at a range of concentrations to
determine the concentration required to achieve 50 % displacement of a reference peptide.
The score is often reported as a 50 % inhibitory concentration (IC 50 ) in nM or as the IC 50
percentile rank.
- MHC **stability** algorithms (e.g. NetMHCstab) predict a related but distinct aspect of the
peptide-MHC complex. While binding prediction models the likelihood that a given peptide
will bind to a specific MHC molecule, stability prediction models how long a peptide-MHC
complex remains stable before dissociating.
- Peptide **cleavage** algorithms (e.g. NetChop) attempt to predict cleavage sites of the
proteasome. These algorithms may be trained on a combination of in vitro protein
degradation and observed MHC ligands. One goal for the interpretation of these algorithms
is to address whether a long peptide, upon cleavage at predicted sites, will retain a core
binding sequence with favorable binding/stability characteristics.
- Peptide **processing** algorithms (e.g. NetMHCpan_EL, MHCflurryEL, BigMHC_EL) are
trained on peptide elution mass spectrometry data and attempt to predict, in aggregate,
peptides that are bound, stable, cleaved and transported appropriately. While conceptually a
superior approach compared to binding affinity predictions, training data remains limited for
some/many HLA alleles.
- **Immunogenicity** prediction algorithms are a rapidly evolving area of neoantigen
prioritization but are significantly limited by availability of training data. The term
immunogenicity is used in many ways, and historically has been used in the context of all of
the above types of prediction algorithms. However, it is most precisely described as relating
to either the prediction of whether a pMHC complex will be recognizable by a TCR or more
broadly whether it is capable of inducing a T cell response in a functioning host immune 
system. Some immunogenicity algorithms directly model the characteristics of the TCR
facing residues of a peptide and these must therefore be interpreted in the context of binding
and presentation predictions. Other immunogenicity algorithms attempt to model all the
requirements of recognition in aggregate. Limited immunogenic true positive/negative
training data makes this last approach challenging at present. Binding and presentation may
be thought of as necessary but not sufficient for immunogenicity.

**Binding Affinity Algorithm Agreement (binding and percentiles scores)**

High agreement in binding affinity score and percentile ranks between multiple algorithms
trained on distinct data with alternative algorithmic approaches may be interpreted as providing
higher confidence in a candidate (with some caveats). In general, the highest agreement
between algorithms is correlated with the amount of training data available for a particular HLA
allele. By contrast, low agreement between algorithms reduces confidence. In cases of low
agreement, particular algorithms may be given greater weight based on published
benchmarking exercises. Furthermore, in the case of poor agreement between binding
algorithms, greater consideration may be given to processing (elution) algorithm scores, based
on the assumption that the independent peptide elution mass spectrometry training data may
provide insight.

Example with high inter-algorithm agreement (binding affinities all < 100 nm)
<img width="1434" alt=" high inter-algorithm agreement" src="https://github.com/user-attachments/assets/d628c5df-c5c8-4b43-8a1d-b9ec684746ec" />

Example with low inter-algorithm agreement (binding affinities range from ~ 300 - 10 , 000 nm)
<img width="1447" alt="low inter-algorithm agreement" src="https://github.com/user-attachments/assets/8d8e85b5-c0b5-42aa-aaea-cab4b0b4d7ba" />

**Processing (aka elution) Algorithm Scores**

Processing (aka elution) algorithm scores are discussed briefly above and generally reported on
a 0 - 1 scale, with higher numbers indicating a higher probability of being presented. The shape
of the distribution of these scores varies from algorithm to algorithm so where possible the


percentile rank should be considered (e.g. considering a percentile score < 2 % to be
acceptable).

**Immunogenicity Algorithm Scores**

Immunogenicity algorithm scores are discussed briefly above and generally reported on a 0 - 1
scale, with higher numbers indicating a higher probability of being recognized by a T cell.
pVACtools currently supports DeepImmuno and BigMHC_IM. DeepImmuno was trained on
immune validation data from IEDB, but not tumor specific immunogenicity data (citation) and
predictions are available for only lengths 9 - 10 and a relatively narrow set of HLA alleles
compared to the binding and presentation algorithms. In general, Immunogenicity prediction
algorithms are considered to be experimental at this stage.

**Cancer Gene Status**
Cancer genes are often defined as those obtained from the Cancer Gene Census. However,
additional review of the known relevance of the specific cancer gene to the specific type of
tumor should also be considered. If a neoantigen corresponding to a cancer gene is not ranked
highly, it may be given additional consideration to see if it can be accepted (e.g. relying on
scores from individual binding, presentation and immunogenicity algorithms instead of median
scores). Cancer Gene neoantigens should be highlighted (e.g. with green background) in the
candidate neoantigen sheet and peptide order sheet.

**Frameshift mutations**

A single frameshift mutation has the potential to give rise to many neoantigens, depending on
how quickly an alternative stop codon occurs in the alternate frame. Frameshift mutations can
introduce considerable complexity to review of neoantigen candidates and selection of long
peptide sequences for the vaccine design. These considerations include:

- Increased possibility of being a false positive somatic variant. Insertion and deletion somatic
variants are perhaps more likely to have quality issues, particularly as the number of
inserted or deleted bases increases (e.g. > 4 ). These variants may be subjected to additional
scrutiny during manual genomic review of the supporting sequence data.
- Miscounting of RNA expression support for the variant allele. Accurately calculating the
variant allele fraction (VAF) for insertion and deletion somatic variants becomes more
difficult as the number of inserted or deleted bases increases (e.g. > 4 ). Manual genomic
review of the supporting sequence data, including consideration of soft-clipped alignments
that support the variant should be performed and the VAF estimate adjusted if necessary.
- Need to adjust length of frameshift sequence analyzed for neoantigens. For performance
reasons we often run pVACtools with a _- -downstream-sequence-length_ setting of 100 amino
acids. While rare, if a frameshift sequence is longer than this before hitting a stop codon,
some neoantigen sequences might be missed. Frameshifts should be reviewed to determine
if a stop codon is reached within the length specified by _- -downstream-sequence-length_.
- Potential for nonsense mediated decay (NMD). In cases where a frameshift variant is
predicted to trigger NMD (e.g. > 50 nucleotides upstream of the last exon-exon junction) the candidate may be considered lower priority, particularly if the RNA allele expression for the
frameshift variant is low compared to the median VAF of all somatic variants for the tumor.
- Additional complexity of proximal variant interpretation. In-phase proximal germline or
somatic variants that occur with the short/long peptide window for a neoantigen candidate
may not be corrected by pVACtools. Proposed short and long peptide sequences should be
investigated to confirm that the impact of any proximal variants on amino acid sequence has
been correctly incorporated into the final peptide sequence.
- Frameshifts that extend the open reading frame beyond the wild type stop codon leading to
a novel C-terminal tail. These candidates may be given high priority because these have
potential for stable expression, avoid NMD concerns, and may lead to long sequences of
tumor specific amino acids.
- Potential to select multiple distinct sequences for the vaccine design. When a frameshift
variant gives rise to a long novel sequence (e.g. > 20 amino acids), it should be examined to
determine if two or more distinct long peptide sequences should be targeted in the final
design. Additional stretches of the frameshift sequence should be examined for qualifying
class I/II targets. Multiple distinct entries may be created in the peptide order form to reflect
distinct synthesis attempts needed to incorporate each into the final design.

Example of a long frameshift peptide sequence (qualifying class I peptide regions in red):
<img width="928" alt="frameshift" src="https://github.com/user-attachments/assets/60b0e1d7-a84a-48b1-b48c-c7e81d57cc81" />


**Dinucleotide mutations**

When dinucleotide variants (DNVs) occur (e.g. chr 14 - 60724007 - 60724008 - GG-TT) pVACtools
will automatically incorporate both nucleotide changes into the predicted short and long peptide
sequences. Neoantigen candidates from these variants can be interpreted similarly to a single
nucleotide variant with one substantial caveat. The VAF for these variants may not be
calculated correctly for these variants in the DNA, RNA or both. Genomics review of the variant
in IGV can be used to determine the correct VAF.

**Anchor Heatmap examination**

The anchor heatmap visualization provided in pVACview reflects an analysis approach
previously described and published ([Xia, et al. Science Immunology, 2023](https://www.science.org/doi/10.1126/sciimmunol.abg2200). A systematic
computational analysis was performed to determine for each HLA allele, what positions of that
allele act to anchor the peptide to the MHC binding groove. Conventional rule of thumb
approaches may consider the first and last few amino acids to be anchoring and the “center” of
the peptide to be TCR facing. A more specific rule of thumb may consider position 2 and the last
position to be the anchor. The anchor heatmap reflects a more systematic approach where the
anchor positions were determined empirically and automatically selected based on a cutoff. By
default the AA positions that contribute 80 % of the anchor potential are defined as anchor
positions. In practice this means that for each allele 1 - 4 positions are generally defined as
anchor positions.

The anchor heatmap (see example below) displays each candidate peptide in the context of the
predicted anchor sites for each HLA allele. The degree to which each position acts as an anchor
is displayed as a normalized anchor score with darker blue highlighting indicating positions that
contribute the most to anchoring. The mutant position is highlighted in red.

If a candidate has a mutation determined to be at an anchor position, the aggretopicity (aka
differential aggretopicity index; DAI) of the candidate should be considered. If the mutation
position corresponds to an anchor, and the binding/presentation of the wild-type peptide is
strong, the candidate should be Rejected. If the mutation position corresponds to an anchor
and the binding/presentation of the wild-type is weak the candidate may be Accepted (assuming
all other criteria are met). If the mutant is not at an anchor position, the wild type
binding/presentation may be ignored. Additional extensive descriptions of this assessment can
be found in the pvactools and pvacview documentation.

<img width="846" alt="Anchor heatmap" src="https://github.com/user-attachments/assets/2e442523-6f27-4733-8489-ee8bad31251f" />

## Genomics Review (post ITB meeting)

Once the ITB meeting is complete and the preliminary set of candidates has been nominated, a
genomics review of all candidates must be completed following the procedures and criteria laid
out below.

## Example Location of Files Needed for Genomics Review

**Variants Review Files**

- **variants.final.annotated.tsv** - > used to verify that Immuno variants were also called by
CLE pipeline

**pVACseq Review Files**

- **$sample-name.all_epitopes.aggregated.tsv** (Class I)
- **$sample-name.all_epitopes.aggregated.metrics.json** (Class I)
- **$sample-name.all_epitopes.aggregated.tsv** (Class II)
- Cancer Gene Census List in TSV format downloaded from Cosmic

**IGV Review Files**

- **annotated.expression.vcf.gz** - > Immuno variant positions to load into IGV
- **normal.cram** - > Normal exome DNA alignments to load into IGV
- **tumor.cram** - > Tumor exome DNA alignments to load into IGV
- **rnaseq/alignments/MarkedSorted.bam** - > Tumor RNA-seq alignments to load into IGV
- **Ensembl 105 _GRCh 38 _UcscGenePred_Custom_Coding.ensGene** - > custom transcript
annotation track containing only transcripts acceptable for neoantigen identification

## Genomics Review Checklist/Principles

The following are high level descriptions of the kinds of detailed review that should be
performed following selection of candidates by the ITB

**Basic data QC review.**
Review relatedness, contamination, FDA QC reports, and purity metrics. Note any red flags.

- Summarize **total number of uniquely mapping reads** generated for Tumor/Normal
exome and Tumor RNA-seq (i.e. not counting duplicates that arise from PCR
amplification of the same source DNA fragment).
- “Unique Mapped Reads” from file:
`qc/fda_metrics/aligned_$sample_dna/table_metrics/$sample_dna_aligned_metrics.txt`
    - Qualitative description of these counts:
        - < 40,000,000 is poor
        - 40,000,000 - 50,000,000 is acceptable
        - 5,000,000 - 100,000,000 is good
        - \> 100,000,000 is excellent
- Summarize **duplication** rates for tumor/normal DNA samples
    - “Mapped Read Duplication” from file:
`qc/fda_metrics/aligned_normal_dna/table_metrics/normal_dna_aligned_metrics.txt`
    - Qualitative description of these rates:
        - \> 75 % is very poor
        - 50 - 75 % is poor
        - 30 - 50 % is acceptable
        - 20 - 30 % good
        - < 20 % is excellent
- Check Somalier results for sample tumor/normal **sample relatedness**.
    - “Relatedness” column from file: “qc/concordance/concordance.somalier.pairs.tsv”
    - Qualitative description of these rates:
        - \> 97.5 % is excellent
        - 95 %- 97.5 % is good
        - 90 - 95 % is concerning
        - < 90 % is very concerning
- Check VerifyBamID results for **contamination** of both tumor and normal samples
    - “FREEMIX” column from file (column 7 ): “qc/$sample/$sample.VerifyBamId.selfSM”
    - Qualitative description of these rates:
        - < 0.025 is good
        - 0.025 - 0.05 is concerning
        - \> 0.05 is very concerning
- Check RNA-seq metrics for **% reads aligning to transcripts**
    - Find the sum of “PCT_CODING_BASES” and “PCT_UTR_BASES” from file:
    `qc/tumor_rna/rna_metrics.txt`
    - Helpful command: “cut -f 17,18 rna_metrics.txt”
        - \> 90 % excellent
        - 75 - 90 % good
        - 50 - 75 % acceptable 
        - < 50 % sub-optimal
- Check for evidence of **end bias** in the aligned RNA seq data using file:
`qc/tumor_rna/rna_metrics.pdf`
    - Visually inspect this plot. If the RNA-seq data is of good quality from intact RNA the
    plot should have a “horse back” shape, representing lower coverage at the beginning
    and end of transcripts but quickly rising and remaining relatively high over the
    majority of the transcript positions. RNA-seq data made from highly degraded RNA
    combined with polyA selection or oligo-dT cDNA priming can have a heavily biased
    distribution instead. Such data can still produce gene expression estimates but may
    be unable to effectively verify expression of some somatic variant alleles.
- Check that the correct RNA strand setting was used in the pipeline YAML file. For
example, the
    - Detected strand file: qc/tumor_rna/trimmed_read_ 1 strandness_check.txt
    - YAML file: `workflow_artifacts/$case_immuno_cloud-WDL.yaml` (grep for “strand”)
    - See further help on the [Troubleshooting README](https://github.com/griffithlab/immuno_gcp_wdl_manuscript/blob/main/Troubleshooting_README.md)
- Summarize total variants called and total neoantigen variants called (i.e. the subset
of variants that could lead to neoantigens). Briefly, total variants is how many somatic
mutations were detected in the tumor. Total neoantigen variants is the subset of these
that lead to possible neoantigen candidates (i.e. those that cause protein coding
changes in known genes). Note that the number of neoantigen candidates selected by
the Immunotherapy Tumor Board will be a small subset of this number.
    - **Total variants** called can be obtained from `variants.final.ansnotated.tsv` (previously
    `pvacseq.annotated.tsv`). Simply the number of rows in this file minus the header.
    - **Total neoantigen variants** called can be obtained from pVACview (total candidates)
    or by looking at the pVACseq aggregated report file to get this number.

To generate a report evaluation of the case's basic QC metrics, run the following commands

```bash
cd $WORKING_BASE

docker run -it --env HOME --env WORKING_BASE -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts:latest /bin/bash

python3 /opt/scripts/get_neoantigen_qc.py -WB $WORKING_BASE -f final_results --yaml $WORKING_BASE/yamls/$CLOUD_YAML
```

**FDA Quality Thresholds**

The following subset of QC values can be obtained from the FDA report tables generated by the
pipeline and used to determine whether the case meets basic data quality criteria as described
in documentation provided to the FDA.

<img width="666" alt="FDA Quality Thresholds" src="https://github.com/user-attachments/assets/af81ee46-3fc1-4c6b-8d7c-c25bda2a71b8" />

To generate a table evaluating the case's FDA metrics, run the following commands
```bash
cd $WORKING_BASE

docker run -it --env HOME --env WORKING_BASE -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts:latest /bin/bash

python3 /opt/scripts/get_FDA_thresholds.py -WB  $WORKING_BASE -f final_results
```

**Tumor type / driver variant review**

Given the reported tumor type do we see evidence for expected somatic driver variants?
Enumerate these in the report and if there are no plausible driver variants, raise this as a
possible red flag. This review can be conducted in pVACview with CancerGeneCensus genes
loaded.

Other strategies for assessing drivers:
- Use the variants TSV file and select all variant types that would directly impact protein
sequence
- Look for variants previous observed in COSMIC or pathogenic according to ClinVar
- Intersect the genes with these mutations with Cancer Gene Census
- Take the candidate genes that result from the previous two queries and search for these
in CBioPortal after selecting studies that match the cancer type of the patient

**HLA allele review for sample/data mixup**
Make sure the HLA alleles as being used in pVACview for the final selection of candidate match
up with those expected for the case based on Optitype/PHLAT predictions from the data and
clinical HLA typing results if available.

**HLA allele review of normal vs tumor**
Check whether the HLA alleles predicted for normal and tumor samples are in agreement. Note
any discrepancies.

To generate a table with the predicted normal and tumor HLA allele, run the following commands
```bash
cd $WORKING_BASE

docker run -it --env HOME --env WORKING_BASE -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts:latest /bin/bash

python3 /opt/scripts/hla_comparison.py -WB $WORKING_BASE
```


**Review multi-algorithm support for strong binding affinity**
We are often using the “lowest” score method to prioritize candidates (the other, more
conservative option would be “median”). This means the peptide with the lowest IC 50 from
*any* algorithm. This leads to more candidates but can lead to the situation where a peptide is
nominated and most algorithms actually predict that it is a poor binder and a minority of
algorithms (perhaps only one) predict that it is a strong binder. In those cases, is the candidate
still acceptable? Is there another peptide where the best score was slightly higher but there is
better agreement across algorithms?
○ MHCnuggets has often been observed to be an outlier. If only this algorithm suggests a
peptide is a strong binder, these candidates should probably be failed in review.
○ NetMHCPan and MHCFlurry have been reported by some benchmarking exercises as
among the better performing algorithms
○ In cases where binding algorithms have high disagreement, additional weight should be
given to presentation algorithms trained on peptide-elution mass spectrometry data

**Compare variant results to CLE pipeline**

**THIS SHOULD BE GENERALIZED OR REMOVED FOR PUBLICATION**

For all variants with neoantigen candidates, check if the variant was also called by the CLE
pipeline. This status is automatically provided in the “VALIDATED” column in the
variants.final.annotated.tsv file produced by the immuno pipeline. The CLE validation status is
also automatically incorporated into the Neoantigen Candidates review spreadsheet.

**IGV Review**
You can SMB mount to obtain access to the case data files on your local computer and create
an IGV session with the following 6 components in order:
1. Final somatic variant VCF from immuno pipeline (annotated.expression.vcf.gz)
2. CLE variant VCF (annotated_filtered.vcf.gz.commented.gz.commented.gz)
3. Normal exome DNA alignments (normal.cram)
4. Tumor exome DNA alignments (tumor.cram)
5. Tumor RNA-seq alignments (rnaseq/alignments/MarkedSorted.bam)
6. Ensembl v 105 transcript annotations (Homo_sapiens.GRCh 38. 105 .sorted.coding.gtf)

Save the session file on storage 1 for others to use. Use this session to perform the following
specific review activities and make note of the findings in the candidate review spreadsheet.

- **Somatic variant review**: For any peptide candidate being considered for inclusion in the
design, manually review the underlying somatic variant in IGV to make sure it appears to
be a real somatic variant. Our [published guidelines](https://www.nature.com/articles/s41436-018-0278-z)
can be used for this step.
- **Proximal variants review**: Check for any nearby variants that might impact prediction of
the peptide sequence. Most of these should be handled automatically by pVACtools if it
was run with the proximal variants option. The definition of “nearby” for the purposes of
proximal variants review is any variant capable of impacting either the predicted class I
epitope sequence ( 8 - 11 amino acids), the class II epitope sequence ( 12 - 18 amino acids)
or the long peptide sequence that will be used for DNA/RNA vector or peptide
sequencing ( 35 - 50 amino acids ). If a variant is NOT near the edge of an exon, the
proximal space to review may be up to 75 bp in either direction (potentially longer for a
frameshift variant). If the variant is near the edge of an exon, the proximal space to
review can be much larger as adjacent exons may be separated by large introns.
Extensive analysis of the importance of proximal variants analysis has been published
([Hundal et al. Nature Genetics, 2019](https://www.nature.com/articles/s41588-018-0283-9)).
- **RNA-seq allele review**: Check for expression of the allele. If the candidate was
selected, but had poor read support for the mutant allele, check whether this could be
explained by a local coverage or alignment issue, where the gene is otherwise well
expressed.
- **Transcript isoform review**: Check if the RNA transcript the candidate peptide is
extracted from is actually expressed. Does the RNA-seq data support the expected
nearby exon-exon junctions? Is there evidence for alternative splicing that might
influence expression of the mutant allele or the structure of the transcript that expresses
it (possibly impacting the presumed peptide sequence)? When performing this review
the following IGV display settings are useful: ( 1 ) “View as pairs” and ( 2 ) “Color
alignments by” - > “first-of-pair strand”. Note that to perform the transcript isoform review,
it will be necessary to load the correct Ensembl GTF of transcripts in IGV. 

**MAKE SURE ENSEMBL GTF IS AVAILABLE**

**Long Peptide Extraction and Annotation**

Using the selected transcript, extract 51-mer peptide sequences that contain the candidate neoantigen. Annotate this sequence with the "best" class-I and class-II binding peptides. The selection of long peptide sequences chosen should reflect the final conclusion of the ITB review and genomics review.

- Use `pvacseq generate_protein_fasta` to create a fasta file that has the long peptide sequences needed.
- For convenience you can create two versions of the Fasta file:
    1. One that has a mutant peptide sequence for every mutation that was evaluated "Accept,Review" during the ITB meeting and based only on the top transcript. This should be close to the final set of sequences you want.
    2. One that has a mutant peptide sequence for every mutation, generated for every transcript. This may be needed to manually adjust the final set of peptides to account for transcript expression, alternative splicing, etc.

An example of how this command would look:
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

## Final Reporting

Generate files for final reporting:

```bash
docker pull griffithlab/neoang_scripts
docker run -it --env WORKING_BASE --env PATIENT_ID -v $HOME/:$HOME/ -v $HOME/.config/gcloud:/root/.config/gcloud griffithlab/neoang_scripts /bin/bash

cd $WORKING_BASE
mkdir manual_review

python3 /opt/scripts/generate_reviews_files.py -reviewed_candidates itb-review-files/*.tsv -peptides generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv  -classI final_results/pVACseq/mhc_i/*.all_epitopes.aggregated.tsv -classII final_results/pVACseq/mhc_ii/*.all_epitopes.aggregated.tsv -samp $PATIENT_ID -o ../manual_review/

python3 /opt/scripts/color_peptides51mer.py -peptides ../manual_review/*Peptides_51-mer.xlsx -samp $PATIENT_ID -probPos C -cIIC50=1000 -cIpercent=2 -cIIIC50=500 -cIIpercent=2 -o ../manual_review/
```

### 1. Peptide Fasta

A plain text fasta file with each peptide sequence. This will be used for final blast checks and can be used with Pepstats to obtain molecular weights for each peptide for the ordering form.

1. Note that if any changes to peptides are made manually during the genomics review (e.g. to select a different reference transcript, or to correct for a complex variant or proximal variant), this file must be updated to have a one-to-one relationship with the Long Peptides Spreadsheet.
2. The starting point for this file is the `annotated_filtered.vcf-pass-51mer.fa` file produced by the `pvacseq generate_protein_fasta` tool described above.
   
   Example location: `$case/generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa`

### 2. Long Peptides Spreadsheet

An annotated peptide sequence file used to create the final vaccine manufacturing order form. The starting point for this spreadsheet is the `.fa.manufacturability.tsv` produced by the `pvacseq generate_protein_fasta` tool described above. 

Example location: `$case/generate_protein_fasta/candidates/annotated_filtered.vcf-pass-51mer.fa.manufacturability.tsv`

The Long Peptide Spreadsheet should provide the following:

#### a. Peptide highlights
ClassI (red text), ClassII (**bold** text) and mutant peptides (underlined) and if using a long peptide vaccine Cysteines (C) indicated with a larger font. Please note the underline convention for these variant types:

1. SNVs - underline the single mutant amino acid
2. In-frame insertions - underline the inserted amino acids
3. In-frame deletions - underline the one amino acid to the left and one amino acid to the right of the deleted amino acid position
4. Frameshift Insertions or Deletions - underline the entire novel amino acid sequence (i.e. from where the sequence starts to diverge from wild type until the end of the peptide sequence)
5. In-frame RNA fusions - underline one amino acid on the left side of the fusion junction and one amino acid on the right side of the fusion junction
6. Frameshift RNA fusions - underline the entire novel amino acid sequence (i.e. from where the sequence starts to diverge from wild type until the end of the peptide sequence)

#### b. Top class I/II HLA alleles
List the best class I and class II allele separated by a "/". List the class I allele if median affinity < 1000 nm OR percentile < 2%. List the class II allele if median affinity < 500 nm OR percentile < 2%. Note that the class I/II peptide should only be red/bold if it meets these criteria.

#### c. Molecular weight
Get the molecular weight for each long peptide sequence by using the EMBOSS Pepstats program.

### 3. Reviewed Candidates Table (MHC Class I)

A spreadsheet with the reviewed class I candidates from pVACview along with any notes on each candidate. Each row corresponds to a single variant leading to potential neoantigens, with the top candidate provided. Each should include the final evaluation: Reject or Accept. Each should also include any detailed notes from the pVACview and IGV reviews that document the rationale for Accepting or Rejecting the candidate.

The starting point for this spreadsheet is the TSV that is exported from pVACview at the end of the review in the Immunotherapy Tumor Board meeting. 

Example storage location and file naming convention: `$case/itb-review-files/mcdb047.revd.Annotated.Neoantigen_Candidates.xlsx`

### 4. Reviewed Candidates Table (MHC Class II)

A spreadsheet with the best class II candidate peptides. Each row corresponds to a single variant leading to potential neoantigens, with the top candidate provided.

The starting point for this spreadsheet is the class II aggregate report from pVACseq that came from the WDL pipeline (or a manual run of pVACseq if that was required). 

Example storage location: `mcdb046/gcp_immuno/final_results/pVACseq/mhc_ii/mcdb046-tumor-exome.all_epitopes.aggregated.tsv`

### 5. Executive Report (aka Genomics Review Report)

A narrative report with an executive summary of key findings from each of the major genomics review steps above. This should be created as each step is completed. In addition to summarizing each step, note any special situations with individual candidates. For example, if a complex variant was manually corrected, if candidates were altered or added to reflect transcript annotation issues or observed alternative splicing, if special consideration was given to driver variants, candidates selected based on class II alone, reference matches, etc. If any peptides were determined or corrected by manual work or ad hoc analysis, describe those in detail here along with visualizations if helpful.l.


## Advanced Topics

**Neoantigen Candidates from RNA Fusions**
The immunogenomics pipeline performs RNA fusion calling with STAR-Fusion and runs
FusionInspector on the resulting fusion candidates. STAR-Fusion predicted fusions are also
annotated with AGfusion and those annotated as having a coding impact are supplied to
pVACfuse which attempts to identify neoantigens. Candidate RNA fusions should be manually
reviewed using FusionInspector and IGV (e.g. review support from soft clipped reads spanning
exon-exon junctions of the fusion). If the fusion support is acceptable (ideally suggesting robust
expression), examine pVACfuse results to determine if any high quality neoantigens are
predicted.

### Fusion Review
Open the fusion inspector html '/gcp_immuno_workflow/rnaseq/fusioninspector_evidence/finspector.fusion_inspector_web.html', this web page will show possible fusions with evidence. A believable fusion would be one with 
- Junction reads + Spanning read counts > 5 and junction reads >= 1
- The fusion is not a read-through
  - Left Chr and Right Chr are different OR chromosome are the same BUT Left Strand and Right Strand are different OR chromosome and strand are the same BUT ABS(Left Pos - Right Pos) < 1,000,000 OR Fusion GeneA Name OR Fusion GeneB Name matches a known fusion driver gene
- The fusion has large anchor support 

We have created [fusion review scripts](https://github.com/kcotto/fusion_review_initial/tree/main) that pull our fusions from the above-listed criteria. 
**SHOULD WE INCLUDE THE REVIEW SCRIPTS? what do you do next? -- test for binding?? Run pvacfuse??**


#### Example: GMB119
<img width="1706" alt="GMB119 fusion inspector" src="https://github.com/user-attachments/assets/cbb288bd-066c-4640-8363-bb47b9df1f4e" />

The PLCG1::LEUTX fusion candidate mentioned above is predicted to produce two possible transcript fusions:

211.PLCG1_LEUTX.ENST00000244007_ENST00000638280.inframe_fusion.72
210.PLCG1_LEUTX.ENST00000244007_ENST00000396841.frameshift_fusion.72

<img width="1507" alt="GBM119 fusion binding" src="https://github.com/user-attachments/assets/93cb366e-16e4-42db-8ecc-6a7fbf00d3c7" />

The frameshift version is not predicted to lead to a good class I peptide, but the inframe version is predicted to give rise to the following peptide with a percentile rank of 2.3% (for HLA-C*08:02).  IGV review of the fusion breakpoints suggests strong support (from soft-clipped reads) on both sides of the fusion in all three tumor regions samples. 

Best peptide: GADKIEGAK 

The fusion CDS sequence from FusionInspector is:
MAGAASPCANGCGPGAPSDAEVLHLCRSLEVGTVMTLFYSKKSQRPERKTFQVKLETRQITWSRGADKIEGAKGPRRYRRPRTRFLSKQLTALRELLEKTMHPSLATMGKLASKLQLDLSVVKIWFKNQRAKWKRQQRQQMQTRPSLGPANQTTSVKKEETPSAITTANIRPVSPGISDANDHDLREPSGIKNPGGASASARVSSWDSQSYDIEQICLGASNPPWASTLFEIDEFVKIYDLPGEDDTSSLNQYLFPVCLEYDQLQSSV*

While the fusion is certainly compelling, the best peptide spans the fusion breakpoint by only a single amino acid. “GADKIEGA” match the wild type PLCG1 sequence. The first amino acid on the LEUTX side, an “E” becomes “K” in the fusion product.  The rest of the sequence from that point is wildtype LEUTX sequence. In other words the amino acid sequence largely presented to the TCR is mostly just wild type PLCG1 sequence.  

MT peptide: GADKIEGAK vs. WT peptide: GADKIEGAE

The mutant and wild type versions of this peptide are not predicted to differ significantly in their binding/presentation.  For this reason, this fusion sequence is not a good candidate. 

#####  Further examples:
https://pvactools.readthedocs.io/en/latest/pvacview/pvacseq_module/pvacseq_vignette.html



**Neoantigen Candidates from Tumor Specific RNA Splicing Events**
A conservative approach to identifying splice neoantigens is to investigate splicing variants that
impact canonical splice acceptor or donor sites (cis regulatory splice variants). In order to
consider a candidate of this type the following criteria should be considered.

- A somatic variant is annotated by VEP as a splice variant (or splice region variant)
- The somatic variant is of high quality (passes manual review as described above)
- The tumor RNA-seq data shows clear evidence of a novel splicing pattern (exon skipping
or alternative acceptor/donor site usage) that makes sense given the splice variant. For
example, a splice acceptor variant is identified and the RNAseq data suggests that the
impacted exon is skipped or an alternative upstream or downstream acceptor site variant
is used. Note that a single somatic splice variant can lead to more than one tumor
specific alternative splicing event.

Assuming the candidate looks good according to these criteria, the protein coding impact of the
alternative splicing can be considered. An alternative transcript can be constructed manually,
translated, and the “neoORF'' can be examined to find novel peptide sequences that can be
processed with pVACbind. Depending on the impact of the splice event the neoantigens may
correspond to (A) a junction peptide sequence arising from exon skipping and inframe
connection of two known exons, (B) a junction peptide arising from alternative acceptor/donor
usage that shortens the transcript and leads to an inframe connection, (C) new inframe peptide
sequence arising from partial intron inclusion, (D) peptide arising from a frameshift induced by
the splicing event.

Note that examination of these impacts and prediction of the neoantigen candidate peptides that
may arise has been implemented as pVACsplice (currently in beta).

**Peptide Modifications to Improve Solubility or other Manufacturing Considerations**
It is possible that a peptide manufacturer may need to add non-wildtype amino acids to the N-
and/or C- terminal end of synthetic long peptides to improve their synthesis, stability or solubility
characteristics. In these instances the modified peptide sequences should be examined with
pVACbind to test whether new strong binding peptides to the patient’s HLA alleles are
inadvertently created and should be excluded.

**Vector Designs**
The purpose of this step is to prepare a construct that contains all the neoantigen peptide
candidates concatenated together with spacers as preferred by a particular manufacturer or as
needed to avoid “junctional epitopes”. These are novel peptide sequences created by the joining
of adjacent candidates. For safety reasons these are screened and any such sequences with
strong MHC binding predictions are avoided. To create a vector design without any such
sequences we use pVACvector which will iterate over possible alternative orderings of the
peptides and test each combination for junctional epitopes to ensure there are none with a
predicted binding affinity less than an accepted cutoff (default 500 nm). If no solution can be
found, the use of spacers (additional non reference amino acid sequences) is attempted and the
iterative process of searching for an optimal ordering is repeated.

The output of pVACvector can be used as input for a DNA, mRNA or adenovirus vaccine
delivery platform.




# Things to ADD
## Add something about running pVACvector
## Fusion Review
## splce mutations?