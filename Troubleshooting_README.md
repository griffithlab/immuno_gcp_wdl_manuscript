# Troubleshooting README

 A pipeline run can take hours to days, and when errors occur, it is not always obvious what is going wrong. Here we will provide a guide on how we would go about troubleshooting a pipeline run and give examples and tips on common errors that we have encountered. 

 # Pipeline Execution Errors

>[!IMPORTANT]  
>Do not delete any Cromwell-executions folders until your pipeline run has excuted sucessesfully. If you deleted these folders, call-caching cannot be utilized and resources (time/money) will be wasted re-executed steps that have already been completed sucessfully.

If you follow these steps and still cannot figure out the error. We encourage posting an issue on GitHub.

## YAML Issues

The most common problem that you can encounter is an error in the YAML file you created. It is important to verify that all the imformation in the YAML file is correct. Tiny errors are easily missed! We have created a YAML checker script to help catch these small, but consequential errors. 


The things the YAML checker script verifies:

- That the sequencing file paths lead to a real file and are not duplicated in the YAML file.
- The readgroup (RG) field is unique between the Normal/Tumor DNA and RNA
- There are no formatting issues, such as improperly commented lines
- If there were problematic amino acids or patient-specific HLA alleles specified
- The stand setting is specified correctly

Run the YAML checker script as follows (we suggest doing this before starting the pipeline run):
   

```bash

docker run -it --env HOME --env GCS_CASE_NAME -v $HOME/:$HOME/ -v /shared/:/shared/ -v $HOME/.config/gcloud:/root/.config/gcloud mgibio/cloudize-workflow:latest /bin/bash

cd $HOME/yamls

python3 /opt/scripts/validate_immuno_yaml.py ${GCS_CASE_NAME}_immuno_cloud-WDL.yaml

```

## Interpreting the Cromwell Log

If the YAML file has been double-checked and confirmed correct, then it is time to further investigate the Cromwell logging file. 

To view the entire cromwell file use this command:

```bash
journalctl -u cromwell
```

### Locate Error

Start at the bottom of the file, where there should be a line that looks like:

```
Aug 18 17:57:49 [VM NAME] java[1710]: 2025-08-18 17:57:49 cromwell-system-akka.dispatchers.engine-dispatcher-7125 INFO  - WorkflowManagerActor: Workflow actor for [WORKFLOW ID] completed with status 'Failed'. The workflow will be removed from the workflow store.
```

The error that caused the pipeline to fail is typically within the lines above this. Depending on the tool, the error message display will look different. Note that if you search the file for key words like 'Error' or 'Fail', you will see a lot of messages containing phrases/command arguments which have to do with error handling. These are not pipeline errors!

Once you find where the error is located, there is usually some indication of which step the error occured at. Sometimes the error message printed in the Cromwell log file is useful, but often it is neccassary to find the `stderr` file for the task from which the error is coming from. Usually there is a line before the error that looks like this: 


```
Aug 18 17:57:45 [VM NAME] java[1710]: Check the content of stderr for potential additional information: gs://[GCS BUCKET]/cromwell-executions/immuno/[WORKFLOW ID]/call-somaticExome/somaticExome/[UNIQUE SUB WORKFLOW ID]/call-detectVariants/detectVariants/[UNIQUE SUB WORKFLOW ID 2]/call-filterVcf/filterVcf/[UNIQUE SUB WORKFLOW ID 23]/call-filterVcfDepth/attempt-3/stderr.
```

### Investigate the stderr file located in the Google Bucket

The error files for all individual tasks are located in the Google Bucket you designated for your results to go. Simply follow the path located in the error file to the stderr file and view it on a web browser. Oftentimes the errors located in this file are much more interpretable and will lead to a conclusion about what the problem is. 

However, if the error is still not interpretable, we begin to do a traceback of pipeline tasks to identify where exactly the problem occurred. In the same folder as the stderr file, explore the other files in the folder, like the `stdout`, `log`, or `script`. The `script` file will let you know what command was run and the inputs needed for that command. Make sure the inputs look correct. **Some common things we have seen are the input files are empty or incredible small, indicating that a process was not completed correctly or was interrupted.**

Note: If `call-caching` was used becasue you resubmitted a workflow on the same VM, there might be a `place_holder.txt` that indicates where the data is being pulled (cached) from. The file will look like this:

```
This directory does not contain any output files because this job matched an identical job that was previously run, thus it was a cache-hit.

Cromwell is configured to not copy outputs during call caching. To change this, edit the filesystems.gcs.caching.duplication-strategy field in your backend configuration.

The original outputs can be found at this location: gs://[GCS BUCKET]/cromwell-executions/immuno/[WORKFLOW ID]/call-rna/rnaseqStarFusion//[UNIQUE SUB WORKFLOW ID]/call-kallisto
```

### Running a step independently 

It might be useful to run a step separately from the rest of the pipeline to try and isolate the issue further. Using the `script` file and the docker from the WDL where that step’s instructions are specfied, you can try and manually rerun the step on VM to try and understand what is failing. 

For example, if we were getting an error during the opitype step, we could find the `script` file used for opitype:

```
gs://[GCS BUCKET]/cromwell-executions/immuno/[WORKFLOW ID]/call-optitype/script
```

Where we can see the command to execute Opitype and the files needed:

```
/bin/bash /usr/bin/optitype_script_wdl.sh /tmp . \

optitype_tumor /mnt/disks/cromwell_root/[GCS BUCKET]/cromwell-executions/immuno/[WORKFLOW ID]/call-somaticExome/somaticExome/[UNIQUE SUB WORKFLOW ID 1]/call-tumorIndexCram/tumor.cram /mnt/disks/cromwell_root/[GCS REFERENCE FILES BUCKET]/human_GRCh38_ens105/aligner_indices/bwamem2_2.2.1/all_sequences.fa 8 50
```

Find the WDL where the opitype execution command is located. It is often helpful to look for tasks on the [analysis-wdls github page](https://github.com/wustl-oncology/analysis-wdls/tree/c0edf02bf3f84766f13df288f3546e499c4bef54) and use search functions to parse through the many WDLs which make up the pipeline. All the specifications for running individual tools are located within the definitions/tools/ directory. So opitype is located at `definitions/tools/optitype_dna.wdl`. 

In this WDL we see the docker image ued for this task: `mgibio/immuno_tools-cwl:1.0.2`

Now we have everything we need to rerun the opitype command, independent of the rest of the pipeline. So on the Google VM:


```
# Copy the files that are needed to rerun the tool

gsutil cp gs://[GCS BUCKET]/cromwell-executions/immuno/[WORKFLOW ID]/call-somaticExome/somaticExome/[UNIQUE SUB WORKFLOW ID 1]/call-tumorIndexCram/tumor.cram .

gsutil cp gs://[GCS REFERENCE FILES BUCKET]/human_GRCh38_ens105/aligner_indices/bwamem2_2.2.1/all_sequences.fa .

# Enter an interactive Docker

docker run -it --env HOME --env GCS_CASE_NAME -v $HOME/:$HOME/ -v /shared/:/shared/ -v $HOME/.config/gcloud:/root/.config/gcloud mgibio/immuno_tools-cwl:1.0.2 /bin/bash

# Run the command and watch for errors

/usr/bin/optitype_script_wdl.sh /tmp . \

optitype_tumor tumor.cram all_sequences.fa 8 50
```

> [!NOTE]  
> Docker will have to be installed on the VM first.

## Running out of Space

Sometimes we encounter errors that are due to the task running out of space during execution. We could see this when sequencing files are larger than expected or if, for example, there are a high amount of mutations called.

Here is an example of what a memory error might look like in the Cromwell logs:

```
Feb 14 06:32:20 [VM NAME] java[1760]: 2025-02-14 06:32:20,715 cromwell-system-akka.dispatchers.engine-dispatcher-91886 

INFO  - WorkflowManagerActor: Workflow [WORKFLOW ID] failed (during ExecutingWorkflowState): java.lang.Except

ion: Task rnaseqStarFusion.kallisto:NA:4 failed. The job was stopped before the command finished. PAPI error code 9. Please check th

e log file for more details: gs://[GCS BUCKET]/cromwell-executions/immuno/[WORKFLOW ID]/call-rna/rnaseqStarFusion/[UNIQUE SUB WORKFLOW ID]/call-kallisto/attempt-4/kallisto.log.
```

As you can see, the error is occruing at kallisto and we are directed to the Kallisto log. If we look at that log, we see that the task tried to pull a certain file onto the task's VM and failed. 

```
Caught ResumableDownloadException (Transfer failed after 23 retries. Final exception: b'[Errno 28] No space left on device') for download of /cromwell_root/[GCS BUCKET]/cromwell-executions/immuno/[WORKFLOW ID]/call-rna/rnaseqStarFusion/[UNIQUE SUB WORKFLOW ID]/call-sequenceToTrimmedFastq/shard-3/sequenceToTrimmedFastq/[UNIQUE SUB WORKFLOW ID 2]/call-trimFastq/glob-[UNIQUE SUB WORKFLOW ID 3]/trimmed_read_2.fastq component 7.

ResumableDownloadException: Transfer failed after 23 retries. Final exception: b'[Errno 28] No space left on device'

CommandException: Some components of /cromwell_root/[GCS BUCKET]/cromwell-executions/immuno/[WORKFLOW ID]/call-rna/rnaseqStarFusion/[UNIQUE SUB WORKFLOW ID]call-sequenceToTrimmedFastq/shard-3/sequenceToTrimmedFastq/[UNIQUE SUB WORKFLOW ID 2]/call-trimFastq/glob-[UNIQUE SUB WORKFLOW ID 3]/trimmed_read_2.fastq were not downloaded successfully. Please retry this download.
```

The most obvious indication that we have a space issue is: ` No space left on device'`. 

To add more space for Kallisto we will edit the WDL for Kallisto on our VM. It may be helpful to go to the [github for the WDLs](https://github.com/wustl-oncology/analysis-wdls/tree/main) to help you find exactly where that task is located within the analysis-wdls file structure. Once you find the task you want to edit, we will change the amount of space needed. Within the WDLs `runtime` section, we passed the space needed to the `disks` attribute. These attributes tell the execution engine (Cromwell) what computational resources your task needs to run. In many of the pipeline's WDLs the space needed is calculated in terms of the size of the files. Ususally this is sufficient, but if we have a space problem, we add just a little extra.

```
# open the WDL you want to edit in a text editor:

vim /shared/definitions/tools/kallisto.wdl

# edit the appropriate line

# for this example I changed line 11 form
# Int space_needed_gb = 10 + round(size(flatten(fastqs), "GB") + size(kallisto_index, "GB"))
# to
# Int space_needed_gb = 30 + round(size(flatten(fastqs), "GB") + size(kallisto_index, "GB"))

# Save your changes

# Refresh the zipped WDLs
source /shared/helpers.sh
refresh_zip_deps

# Resubmit workflow
submit_workflow /shared/analysis-wdls/definitions/$WORKFLOW $YAML 

```

> [!IMPORTANT]
> Make sure to run `refresh_zip_deps` which will generate a zipped version of the WDLs that will be passed to Cromwell to execute. If you don't, the changes that you made will not be applied.

> [!TIP]
> Space vs Memory: When looking at the WDLS runtime attributes, you might be considering whether you need to fix the amount of disks or the amount of memory. Keep in mind that disks refer to the space you have to store your files, while memory refers to the amount "space" you need to do your computational task. Memory will rarely have to be changed.

# Errors in Results

Sometimes the problem is not that the pipeline does not complete sucessfully but that the results are not as expected. Here we will go through some possible scenarios on incorrect results and what to do about it.

## Adjusting Other Parameters in the YAML

There are many thresholds set in the YAML file that have been tested over hundreds of samples, but sometimes adjusting them is necessary in most cases becasue there is sparse data. For example, you have pancreatic tumor data and expect to see a KRAS mutation, but the tumor content is low.

To reduce variant filtering, we can change the VAF cutoff for Varscan, the VAF cutoff for the fp filter, and also the log likelihood ratio test threshold to allow lower VAF variants through. 

```
#reduce filtering stringency
immuno.varscan_min_var_freq: 0.025   # default is 0.05
immuno.fp_min_var_freq: 0.025   # default is 0.05
immuno.filter_somatic_llr_threshold: 3.5 # default is 5
```

There may be cases where the peptide candidates are sparse, and thresholds need to be adjusted for to let more candidates through. For example:

```
immuno.aggregate_inclusion_binding_threshold: 2500 # set between 1500 and 5000 depending on how many candidates are expected
immuno.downstream_sequence_length: 100 # may need to be increased if a very long frameshift is seen
immuno.min_ffpm_level: 0.025 # use a level as low as this for sensitive detection of fusions
```

## The strand settings

A typical thing that the qc review can reveal is that the wrong strand setting was used in the YAML file. It is important to know the correct strand setting to determine the expression levels of a gene, especially in cases of overlapping genes. So with an unstranded protocol, we do not map RNA reads to DNA by strand information. When the strand protocol is unknown, setting the strand setting in the yaml to unstranded is safest because it essentially does not consider the directionality of the reads during mapping. However, this does result in a loss of information. If you set your yaml to first strand when the protocol is second or unstranded, you could potentially lose whole genes depending on the overlap.

The pipeline's inferred strandedness result can be viewed after you pull results (`./final_results/qc/tumor_rna/trimmed_read_1strandness_check.txt`). This file contains this at the bottom:

```
Fraction of reads failed to determine: 0.0761
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0380
Fraction of reads explained by "1+-,1-+,2++,2--": 0.8860
Over 75% of reads explained by "1+-,1-+,2++,2--"
Data is likely RF/fr-firststrand
```

This says that ~4% of reads are explained by the second strand (`1++,1--,2+-,2-+`) and ~89% of reads are on the first strand (`1+-,1-+,2++,2--`).

The first-strand method sequences the first synthesized strand of cDNA, which is complementary to the original RNA.

Read 1 should be on the OPPOSITE strand of the gene.
1+-: Read 1 is forward (+), gene is reverse (-).
1-+: Read 1 is reverse (-), gene is forward (+).
Read 2 should be on the SAME strand as the gene.
2++: Read 2 is forward (+), gene is forward (+).
2--: Read 2 is reverse (-), gene is reverse (-).

If the sequencing provider does not make it clear what strand your data is, it is okay to rely on this result. We recommend rerunning the pipeline with the correct strand to ensure optimal results.

## End Bias

Visually inspect the plot located at `qc/tumor_rna/rna_metrics.pdf`. If the RNA-seq data is of good quality from intact RNA, the plot should have a “horseback” shape, representing lower coverage at the beginning and end of transcripts but quickly rising and remaining relatively high over the majority of the transcript positions. RNA-seq data made from highly degraded RNA, combined with polyA selection or oligo-dT cDNA priming, can have a heavily biased distribution instead. Such data can still produce gene expression estimates, but will be unable to effectively verify the expression of somatic variant alleles.

An example of a good end bias plot:

![5120-30 end-bias horseshoe GOOD](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/2cea9675-5f30-45a2-b6ad-8afe458c7099)


An example of a poor end-bias plot. This case was run with the YAML set to "unstranded" when the detected strandedness of the RNA data was unstranded. **is this the reson of this plot???**

![jlf-100-067 end-bias BAD](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/649f9477-2d85-450c-8257-282c3be60abf)

## Problems Encountered Upon Inspection in IGV 

To fully verify neoantigen peptide candidates, visual inspection in IGV is crucial. Often, IGV review reveals sequencing read issues or algorithmic mistakes that can only be caught visually. Here are a few examples of what you may find.

### High Duplicaton rate

High duplication rate appeared to create situations where all variant support came from a set of identical/duplicate read alignments. In this screenshot, we see a series of reads that look identical in length, which makes them appear to be sequencing artifacts.

![SLF2 variant high duplication rate](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/e5baa217-3cb0-4b4f-9ea4-f21dedee0133)

### Multiple RNA Splicing

Sometimes a candidate has high read support, but when inspecting the RNA read splicing (using a sashmi plot in IGV), we see that the RNA supports splicing to different transcripts. Sometimes the peptide candidate is on a transript that is not as supported as another, or the splicing pattern suggest that the sequence would give rise to a different peptide.

> [!NOTE]  
> **maybe a better picture with the transcripts**

![CPEB2 alternative splicing](https://github.com/evelyn-schmidt/immuno_gcp_wdl_manuscript/assets/57552529/d341bc01-f80c-4817-881f-1695380f4085)

### Alterations to the best peptide sequences

> [!NOTE]  
> Unsure about this section

Visual inspection often reveals that the best peptide, or more commonly 51mer sequence generated for pepetide manufactoring, need to be changed to fit what the data shows. For example, there might be a heterozygous germline variant that changes the sequence but was not caught by the pipeline's germline variant detector.

Changes to the best peptide need to be checked by rerunning pVACseq to make sure that the peptide is still a strong and unproblematic binder. See the pVACseq documentation [here](https://pvactools.readthedocs.io/en/latest/pvacseq.html).

## Google Cloud Errors

### Bucket permissions

```
Aug 05 15:59:34 sidi-immuno-jlf-100-128 java[1744]: Failed to evaluate input 'data_size' (reason 1 of 1): [Attempted 1 time(s)] - StorageException: cromwell-server@jlf-rcrf.iam.gserviceaccount.com does not have storage.objects.get access to the Google Cloud Storage object. Permission 'storage.objects.get' denied on resource (or it may not exist).
```