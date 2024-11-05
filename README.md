# CUT & TAG pipeline 

This is a Cut & Tag processing pipeline for the R-loop identification project collaboration between the Kinkley lab & Vingron Lab.
The raw data was produced by Serkan Meydaneri. The pipeline's goal is to produce BigWig files and Peak files from the raw
Cut & Tag data.

The pipeline consists of modules (bash files for each task) that are then called by a main script `cut_and_tag_pipeline.sh`.
The modules are: 

1. Quality Control: Apply FastQC to the data `1_QC_Preprocessing.sh`.
2. Alignment: Mapping of CnT reads to the hg38 human genome and spike-in genome (`2_Alignment.sh` and `2bis_SpikeInAlignment.sh`).
3. Assessment of alignment: Scripts to assess number of duplicates, sequencing depth, fragment size etc.
4. Filter, conversion and Reproducibility Assessment: Filtering of low quality reads, convert files to bam bed and bw and assess reproducibility
between replicates.
5. IgG substraction and scaling
6. Peak calling with SEACR

## How to run 

0. Create the project folder which will hold all your processed files and 
download the pipeline from github
1. Create the necessary files
The pipeline needs the following csv files

- experiment_summary.csv

| SampleName | File |
|:------------:|:--------------:|
| H3K4me3_1_1  | H3K4me3_1_R1.fastq.gz | 
| H3K4me3_1_2  | H3K4me3_1_R2.fastq.gz |

- experiment_summary_align_formatted_example.csv

| SampleName | FileNameRep1 | FileNameRep2 |
|:------------:|:--------------:|:--------------:|
| H3K4me3_1  | H3K4me3_1_R1.fastq.gz | H3K4me3_1_R1.fastq.gz |
| H3K4me3_2  | H3K4me3_1_R2.fastq.gz | H3K4me3_1_R2.fastq.gz |

A table with two replicates of H3K4me3. FileNameRep1 and FileNameRep1 are for the paired-end
reads.

- experiment_summary_peaks.csv

| SampleName | Control |
|:------------:|:--------:|
| H3K4me3_1  | IgG_1 |
| H3K4me3_2  | IgG_2 |

Group each experiment with the control experiment you want to use for peak calling
and for IgG substraction (optional).

2. Copy the file cut_and_tag_pipeline_tutorial.sh or cut_and_tag_pipeline_example.sh, 
change its name and the following variables: 
- PROJECTROOT
- RAWROOT
- PIPELINE
- EXPSUMMARY
- EXPSUMMARYAL
- EXPSUMMARYPEAKS
- DUPREMOVE
- IGGNEED

Explanations on the variables are inside both files

3. Run pipeline on the servers

Open a terminal on your computer and type the following:
```
hostconfig --list public #list all of our available servers
ssh myserverofchoice # change to name of server
cd /my/chosen/folder # change to your folder
mv Cutntag pipeline # change the name of the pipeline folder
cd pipeline
bash cut_and_tag_pipeline_example.sh > cut_and_tag_pipeline_example.out 
#> means you send the output of the previous command to the .out file. You can then check the output of the pipeline in this file
#IMPORTANT!! check the .out file for any errors. If any problems arise on your runs this will tell me what happened.
history #shows you the recent commands
```

## Meaning of file suffixes

The pipeline produces many different files in SAM, BAM, BED, BEDGRAPH and BIGWIG 
formats (explanation of formats: https://genome.ucsc.edu/FAQ/FAQformat.html). 
Meanings of the file suffixes: 
- .normalized: Spike-in normalized file
- .rmDup: Deduplicated file
- .substracted.igg: File with substracted igG

## Naming recommendation

We recommend using "_" to specify the replicates (e.g. H3K4me3_1), if you have 
knockdowns or over-expressions please use PROTEINKD_1 or PROTEIN.KD_1 but not
PROTEIN_KD_1.

