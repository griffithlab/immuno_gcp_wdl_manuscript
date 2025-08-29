# Troubleshooting README

Here we will give examples of problems that you may encounter during the ITB and immunogenomics review process


## YAML Issues

The most common problem that you can encounter is an error in the YAML file you created. 

- file paths being wrong
- not matching sample IDS
- strand setting (this won't mess up the pipleine but it will mess up results)

We created a yaml checker script to help catch these small, but consequential errors. 
   
```bash
docker run -it --env HOME --env GCS_CASE_NAME -v $HOME/:$HOME/ -v /shared/:/shared/ -v $HOME/.config/gcloud:/root/.config/gcloud mgibio/cloudize-workflow:latest /bin/bash

cd $HOME/yamls

python3 /opt/scripts/validate_immuno_yaml.py ${GCS_CASE_NAME}_immuno_cloud-WDL.yaml
```

## Interpretting the Cromwell Log

If the YAML file has been double-checked and confirmed correct, then it is time to further investigate the cromwell logging file. 

To view the entire cromwell file use this command:
```bash
journalctl -u cromwell
```

Start at the bottom of the file where there should be a line which looks like:

```
Aug 18 17:57:49 [VM NAME] java[1710]: 2025-08-18 17:57:49 cromwell-system-akka.dispatchers.engine-dispatcher-7125 INFO  - WorkflowManagerActor: Workflow actor for [WORKFLOW ID] completed with status 'Failed'. The workflow will be removed from the workflow store.
```

The error which caused the pipeline to fail is typically within the lines above this. Depending on the tool, the error message display will look different. Note that if you search the file for key words like 'Error' or 'Fail' you will see alot of messages containing phases/command arguments which are to do with error handling. These are not pipeline errors. 

Once you find where the error is located there is usaully some indication of what step the error is coming from. Sometimes the error message printed in the cromwell log file is useful, but often it is neccassary to find the stderr file for the task which the error is coming from. Ususally there is a line before the error that looks like this: 

```Aug 18 17:57:45 [VM NAME] java[1710]: Check the content of stderr for potential additional information: gs://[GCS BUCKET]/cromwell-executions/immuno/[WORKFLOW ID]/call-somaticExome/somaticExome/[UNIQUE SUB WORKFLOW ID]/call-detectVariants/detectVariants/[UNIQUE SUB WORKFLOW ID 2]/call-filterVcf/filterVcf/[UNIQUE SUB WORKFLOW ID 23]/call-filterVcfDepth/attempt-3/stderr.```

The error files for all individual tasks are located in the google bucket you designated your results to go to. Simply follow the path located in the error file to that stderr and view it on the web browser. Often times the errors located in this file are much more interpretable and will lead to a conclusion about what the problem is. 

However, if the error is still not intreptable, we begin to do a traceback of pipeline tasks to identify where exactly the problem occured. In the same folder as the stderr file, explore the other files in the folder like the `stdout`, `log`, or `script`. The `script` file will let you know what command was run and inputs needed for that command. Make sure the inputs look correct. Some common things we have seen are the input files are empty or incredible small, indicating that a process was not completed correctly or was interrupted. 

Note: If the `call-caching` was used, if you resubmitted a workflow on the same VM, there might be a `place_holder.txt` that indicated where the data is being pulled (cahced) from. The file will look like this:

```
This directory does not contain any output files because this job matched an identical job that was previously run, thus it was a cache-hit.
Cromwell is configured to not copy outputs during call caching. To change this, edit the filesystems.gcs.caching.duplication-strategy field in your backend configuration.
The original outputs can be found at this location: gs://jlf-100-037/cromwell-executions/immuno/598fac4f-d678-49f6-a8a7-cbc9cb7adcc9/call-rna/rnaseqStarFusion/07a05ef5-4fd7-4ba2-8474-3189ad246480/call-kallisto
```

It might be useful to run a step separetly from the rest of the pipeline to try and isolate the issue further. Using the `script` file and the docker from the 

### Example: Sample IDs Not Matching
```
Aug 18 17:57:45 eve-immuno-jlf-100-132 java[1710]: 2025-08-18 17:57:45 cromwell-system-akka.dispatchers.engine-dispatcher-7117 INFO  - WorkflowManagerActor: Workflow ba1ded49-c6b4-42cf-8127-aa0d51ff0883 failed (during ExecutingWorkflowState): Job filterVcf.filterVcfDepth:NA:3 exited with return code 1 which has not been declared as a valid return code. See 'continueOnReturnCode' runtime attribute for more details.
Aug 18 17:57:45 eve-immuno-jlf-100-132 java[1710]: Check the content of stderr for potential additional information: gs://jlf-100-132/cromwell-executions/immuno/ba1ded49-c6b4-42cf-8127-aa0d51ff0883/call-somaticExome/somaticExome/b1985ce0-ded1-4376-a6bb-82c2fe159fe0/call-detectVariants/detectVariants/20e4c9ce-4441-4a83-bc71-213f8c40a6be/call-filterVcf/filterVcf/bb128705-b12b-428a-a80e-5ccd48061891/call-filterVcfDepth/attempt-3/stderr.
Aug 18 17:57:45 eve-immuno-jlf-100-132 java[1710]:  [First 3000 bytes]:Traceback (most recent call last):
Aug 18 17:57:45 eve-immuno-jlf-100-132 java[1710]:   File "/usr/bin/depth_filter.py", line 105, in <module>
Aug 18 17:57:45 eve-immuno-jlf-100-132 java[1710]:     main()
Aug 18 17:57:45 eve-immuno-jlf-100-132 java[1710]:   File "/usr/bin/depth_filter.py", line 93, in main
Aug 18 17:57:45 eve-immuno-jlf-100-132 java[1710]:     elif(depth < args.minimum_depth):
Aug 18 17:57:45 eve-immuno-jlf-100-132 java[1710]: TypeError: '<' not supported between instances of 'NoneType' and 'int'
Aug 18 17:57:49 eve-immuno-jlf-100-132 java[1710]: 2025-08-18 17:57:49 cromwell-system-akka.dispatchers.engine-dispatcher-7125 INFO  - WorkflowManagerActor: Workflow actor for ba1ded49-c6b4-42cf-8127-aa0d51ff0883 completed with status 'Failed'. The workflow will be removed from the workflow store.
```
### Space Issue

```
Jan 29 21:36:11 eve-immuno-jlf-100-082 java[12798]: 2025-01-29 21:36:11,916 cromwell-system-akka.dispatchers.engine-dispatcher-22396 INFO  - WorkflowManagerActor: Workflow 9e157c36-5ce7-42f4-8f16-6c1031486475 failed (during ExecutingWorkflowState): java.lang.Exception: Task mutect.mutectTask:24:4 failed. The job was stopped before the command finished. PAPI error code 10. The assigned worker has failed to complete the operation

Jan 29 21:36:11 eve-immuno-jlf-100-082 java[12798]: java.lang.Exception: Task mutect.mutectTask:3:4 failed. Job exit code 247. Check gs://jlf-rcrf-immuno-outputs/cromwell-executions/immuno/9e157c36-5ce7-42f4-8f16-6c1031486475/call-somaticExome/somaticExome/c822eef7-65c8-4660-8d39-a2afe1ebc6d0/call-detectVariants/detectVariants/43b22f5f-7b8f-4250-81ce-c28c806d5a4b/call-mutect/mutect/52f34299-cd9b-427f-97e5-6f997f7ee85c/call-mutectTask/shard-3/attempt-4/stderr for more information. PAPI error code 9. Please check the log file for more details: gs://jlf-rcrf-immuno-outputs/cromwell-executions/immuno/9e157c36-5ce7-42f4-8f16-6c1031486475/call-somaticExome/somaticExome/c822eef7-65c8-4660-8d39-a2afe1ebc6d0/call-detectVariants/detectVariants/43b22f5f-7b8f-4250-81ce-c28c806d5a4b/call-mutect/mutect/52f34299-cd9b-427f-97e5-6f997f7ee85c/call-mutectTask/shard-3/attempt-4/mutectTask-3.log.
```

## Google Cloud Errors

### Bucket permissions

```
Aug 05 15:59:34 sidi-immuno-jlf-100-128 java[1744]: Failed to evaluate input 'data_size' (reason 1 of 1): [Attempted 1 time(s)] - StorageException: cromwell-server@jlf-rcrf.iam.gserviceaccount.com does not have storage.objects.get access to the Google Cloud Storage object. Permission 'storage.objects.get' denied on resource (or it may not exist).
```

### Google Cloud Outage

```
Jun 13 03:10:40 eve-immuno-jlf-100-119 java[1788]: We encountered an internal error. Please try again.
Jun 13 03:10:40 eve-immuno-jlf-100-119 java[1788]: Caused by: java.io.IOException: Could not read from gs://evelyn-jlf-immuno/cromwell-executions/immuno/d05c070f-df
1c-4917-9d52-545b057c53c2/call-rna/rnaseqStarFusion/

```


## Common Problems Encountered During Immunogenomics Review

### The strand settings
A typical thing that the qc review can reveal is that the wrong strand setting was used in the yaml file. It is important to know the correct strand setting to determine the expression levels of gene, especially in cases of overlapping genes. So with an unstranded protocol, we do not map RNA reads to DNA by strand information. When the strand protocol is unknown, setting the strand setting in the yaml to unstranded is safest because it essentially does not consider the directionality of the reads during mapping. However, this does result in a loss of information. If you set your yaml to first strand when the protocol is second or unstranded, you could potentially lose whole genes depending on overlap...

**how to manually infer strandedness**

### End Bias

Visually inspect the plot located at `qc/tumor_rna/rna_metrics.pdf`. If the RNA-seq data is of good quality from intact RNA the plot should have a “horseback” shape, representing lower coverage at the beginning and end of transcripts but quickly rising and remaining relatively high over the majority of the transcript positions. RNA-seq data made from highly degraded RNA combined with polyA selection or oligo-dT cDNA priming can have a heavily biased distribution instead. Such data can still produce gene expression estimates but will be unable to effectively verify the expression of somatic variant alleles.

An example of a good end bias plot:
![5120-30 end-bias horseshoe GOOD](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/2cea9675-5f30-45a2-b6ad-8afe458c7099)

An example of a poor end-bias plot. This case was run with the yaml set to "unstranded" when the detected strandedness of the RNA data was unstranded. **is this the reson of this plot???**
![jlf-100-067 end-bias BAD](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/649f9477-2d85-450c-8257-282c3be60abf)


###  Assessing Read quality

the high duplication rate appeared to create situations where all variant support came from a set of identical/duplicate read alignments
![JLF-100-060 SLF2 variant high duplication rate](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/e5baa217-3cb0-4b4f-9ea4-f21dedee0133)

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

### Alterations to best peptide that needs reanalysis of binding

### Other Examples -- TO REPLACE LATER

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

#### JLF-100-039_mcdb045
The SLC9A6 on the other hand has a downstream missense variant (almost 4 bases away) that alters the expected best peptide, thus we do not have information of binding affinity, etc. for what would be the best new peptide. This candidate could potentially be rescued with additional effort if needed.

#### JLF-100-037_mcdb041-original
Results from pVACfuse were used to include a peptide for a GFPT1::ENOX2 fusion in this tumor. The peptide sequence corresponds to the tumor-specific frameshift sequence created by the fusion event.

The ALK portion of the EML4::ALK fusion was extracted and used with pVACbind to nominate candidates for that driver event. In this case, while the ALK sequence is the portion presumably amplified/activated by the fusion event, the actual ALK sequence is NOT tumor specific.  It is simply the wild-type ALK sequence, from exon 20 to the end of the protein that is fused as the 3’ component of the EML4::ALK fusion. 

We also examined a fusion prediction for ZNF92::TDRD9 which had 44 junctions reads of support and looks real by manual review in fusion inspector and IGV. However this fusion is predicted to join the first exon of ZNF92 onto the second exon of TDRD9 and this is predicted to lead to an almost immediate stop codon.  The overall ORF would only be 10 AA or so and this is unlikely to be translated.

We also examined a fusion prediction for KDM5C::KMT2C which had 34 junction reads and 2 spanning reads of support.  It is predicted to lead to an almost immediate stop codon.  The neoORF segment would be: MLFHGCLS.  This truncation would be a drastic shortening of the normal KMT2C protein.  The only thing that is predicted to be a good binder is: MSHGVPMLF.  In other words only incorporating the first 3 AA of the novel sequencing arising from the fusion.  This is probably not worth targeting.


##### Leidos 5120-29 
SLC1A7 candidate (original status: Review) was dropped during the manual review. Short explanation: the Class I peptide which doesnt have reference match has bad binding affinity, Class II peptide has reference match. Long explanation: This candidate has 2 transcript sets (ENST00000620347.5 and ENST00000371494.9) , both transcript sets have a matched portion (LQALLIVL) with proteome reference. Class II peptide (QALLIVLATSSSSA, from ENST00000371494.9) has a complete overlap with the reference, thus was rejected. Class I peptide for HLA-C*03:04 (VLATSSSSATL, from ENST00000371494.9) can still be used (if we cut the 51 mer to exclude the reference match portion), however this peptide is a bad binder (median IC50 greater than 2,000 nM, with multiple algorithms reports high IC50). Class I peptide for HLA-A (ILQALLIVL from ENST00000620347.5) has good binding affinity but has a strong reference match.  






