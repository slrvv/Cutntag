################################################################################
# cut & tag analysis pipeline R-loop identification project                    #
# 1. Preprocessing & QC                                                        #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: Preprocessing of the cut & tag data FASTQC and possible merging of  #
# replicates.                                                                  #
################################################################################

#----------------------------Paths to project folders--------------------------#
### command line arguments
#rep1: raw fastq file for replicate 1 of the experiment
#rep2: raw fastq file for replicate 2 of the experiment
#projPath: project path folder containing all analysis on this particular
# experiment.
#name: short name abbreviation to identify experiment (i.e. H3K27ac)
rep1=$1
rep2=$2
projPath=$3
name=$4
ref=$5

#-------------------------Alignment Script-------------------------------------#

cores=8

#mkdir -p ${projPath}/alignment/sam/bowtie2_summary
#mkdir -p ${projPath}/alignment/bam
#mkdir -p ${projPath}/alignment/bed
#mkdir -p ${projPath}/alignment/bedgraph

## Build the bowtie2 reference genome index if needed:
## bowtie2-build path/to/hg38/fasta/hg38.fa /path/to/bowtie2Index/hg38

bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 \
	-X 700 -p ${cores} -x ${ref} -1 $rep1 -2 $rep2 \
	-S ${projPath}/alignment/sam/${name}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${name}_bowtie2.txt
#echo "bowtie not working"
