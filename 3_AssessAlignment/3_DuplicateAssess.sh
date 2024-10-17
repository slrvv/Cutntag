################################################################################
# cut & tag analysis pipeline R-loop identification project                    #
# 3. Alignment assessment                                                      #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: Assessing number of duplicates                                      #
################################################################################

#----------------------------Paths to project folders--------------------------#
### command line arguments
#rep1: raw fastq file for replicate 1 of the experiment
#rep2: raw fastq file for replicate 2 of the experiment
#projPath: project path folder containing all analysis on this particular
# experiment.
#name: short name abbreviation to identify experiment (i.e. H3K27ac)
picard=/scratch/ngsvin/bin/Picard2.21.1/picard.jar
picardCMD="java -jar $picard"
projPath=$1
name=$2
ref=$5

#-------------------------------Dup Script-------------------------------------#

#Technically one shouldn't remove duplicates in cutntag data. But we will assess
#whether in this case it makes sense or not.


mkdir -p $projPath/alignment/removeDuplicate/picard_summary

## Sort by name
$picardCMD SortSam I=$projPath/alignment/sam/${name}_bowtie2.sam \
O=$projPath/alignment/sam/${name}_bowtie2.sorted.sam \
SORT_ORDER=queryname

## mark duplicates
$picardCMD MarkDuplicates I=$projPath/alignment/sam/${name}_bowtie2.sorted.sam \
O=$projPath/alignment/removeDuplicate/${name}_bowtie2.sorted.dupMarked.sam \
METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${name}_picard.dupMark.txt

## remove duplicates
$picardCMD MarkDuplicates I=$projPath/alignment/sam/${name}_bowtie2.sorted.sam \
O=$projPath/alignment/removeDuplicate/${name}_bowtie2.sorted.rmDup.sam \
REMOVE_DUPLICATES=true \
METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${name}_picard.rmDup.txt
