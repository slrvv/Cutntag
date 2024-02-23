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


