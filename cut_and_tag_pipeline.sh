################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# all modules together                                                         #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: This script puts all module scripts together (one script per step)  #
# for the analysis of the Cut&tag data. Paired-end sequencing and no Spike-In  #
################################################################################


#--------------------------Project Paths---------------------------------------#
#project folder
PROJECTROOT=/project/ChromGroup/Serkan_Project/cut_and_tag_rloops
#rootfolder to the raw experimental data
RAWROOT=/project/solexawork/pipelined/230920_SN435_20_AACTTJMM5/demux_20230921_SN435-AACTTJMM5_std-8d1/Project_V
PIPELINE=/project/ChromGroup/Serkan_Project/cut_and_tag_rloops/pipeline
#csv file with all of the experiments and a short id name
EXPSUMMARY=$PROJECTROOT/experiment_summary.csv
#csv file with all of the experiments but replicates are given in columns
#more suitable for the alignment step
EXPSUMMARYAL=$PROJECTROOT/experiment_summary_align_formatted.csv
EXPSUMMARYALREDO=$PROJECTROOT/experiment_summary_alignment_redo.csv

#-------------------------Module Paths-----------------------------------------#
#bash file that does the first step QC
QC=$PIPELINE/1_QC_Preprocessing.sh
#bash file for alignment
ALIGN=$PIPELINE/2_Alignment.sh
#bash file for spike in calibration
SPIKEIN=$PIPELINE/2bis_SpikeInAlignment.sh
#makes a summary of seq depth
SUMMARYAL=$PIPELINE/3_MappingSummary.R
#assess duplicates
DUP=$PIPELINE/3_DuplicateAssess.sh
#Fragment assessment
FRAG=$PIPELINE/3_FragAssess.sh
#summary of fragment and duplicates assessment
DUPSUM=$PIPELINE/3_DuplicateAssessSumary.R
FRAGSUM=$PIPELINE/3_FragAssessSumary.R
SPIKEINSUM=$PIPELINE/3_SpikeinSummary.R
#Filterandconvert
FILTERCONV=$PIPELINE/4_FilterAndConvert.sh
#Replicate reproducibility assessment
REPREPRO=$PIPELINE/4_ReplicateReproducibility.R
#-------------------------1. Quality Control-----------------------------------#

#Quality control each experiment. EXPSUMMARY is a table containing the file names
#for each experiment and a summary name. The FASTQC results are saved under
#the FastQC folder with subfolder name after the short experiment summary name
# split[0] is summary name split[1] is file name
# echo "Quality Control with FASTQC"
# echo " "
# while read line ; do
#     set $line
#     IFS=$','; split=($line); unset IFS;
#     if [[ $accountId == "Sample" ]]; then
#         continue
#     fi
#     echo "Name of experiment ${split[0]}"
#     echo "Name of file ${split[1]}"
#     bash ${QC} $PROJECTROOT/FastQCResults/${split[0]} $RAWROOT/${split[1]}
#     echo " "
# done < <(tail -n +2 $EXPSUMMARY)
# 
# #---------------------------2. Alignment---------------------------------------#
# 
# #Alignment for each experiment, each experiment has 2 replicates EXPSUMMARYAL
# #has one row per experiment and columns per replicate. split[0] is the name
# #of the experiment split[1] is first replicate and split[2] is the second.
# #The experiments were repeated by two different people which is marked by the _1
# #or _2 in the name (split[0])
# cores=8
# ref=/project/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index_Large/genome
# echo "Alignment to hg38"
# echo " "
# mkdir -p $PROJECTROOT/alignment/sam/bowtie2_summary
# mkdir -p $PROJECTROOT/alignment/bam
# mkdir -p $PROJECTROOT/alignment/bed
# mkdir -p $PROJECTROOT/alignment/bedgraph
# while read line ; do
#      set $line
#      IFS=$','; split=($line); unset IFS;
#      if [[ $1 == "SampleName" ]]; then
#          continue
#      fi
#      echo "Name of the sample ${split[0]}"
#      echo "Name of R1 ${split[1]}"
#      echo "Name of R2 ${split[2]}"
#      bash $ALIGN $RAWROOT/${split[1]} $RAWROOT/${split[2]} $PROJECTROOT ${split[0]} $ref
#      echo " "
# done < <(tail -n +2 $EXPSUMMARYAL)
# 
# #Redo alignments that failed
# #cores=8
# #ref=/project/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index_Large/genome
# #echo "Alignment to hg38"
# #echo " "
# #mkdir -p $PROJECTROOT/alignment/sam/bowtie2_summary
# #mkdir -p $PROJECTROOT/alignment/bam
# #mkdir -p $PROJECTROOT/alignment/bed
# #mkdir -p $PROJECTROOT/alignment/bedgraph
# #while read line ; do
# #    set $line
# #     IFS=$','; split=($line); unset IFS;
# #     if [[ $1 == "SampleName" ]]; then
# #         continue
# #
# #     fi
# #     echo "Name of the sample ${split[0]}"
# #     echo "Name of rep 1 ${split[1]}"
# #     echo "Name of rep 2 ${split[2]}"
# #     bash $ALIGN $RAWROOT/${split[1]} $RAWROOT/${split[2]} $PROJECTROOT ${split[0]} $ref
# #     echo " "
# #done < <(tail -n +2 $EXPSUMMARYALREDO)
# 
# 
# #-------------------3. Alignment Summary and assessment------------------------#
# 
# #Produce a table summarizing the sequencing depth of all of the alignments
# echo "Make Alignment Summary"
# Rscript $SUMMARYAL $PROJECTROOT
# echo " "
# #Create summary plots for the alignments using the 3_MappingPlots.R script
# 
# #Assess the number of duplicates
# 
# echo "Assess Duplicates"
# while read line ; do
#     set $line
#     IFS=$','; split=($line); unset IFS;
#     if [[ $1 == "SampleName" ]]; then
#         continue
# 
#     fi
#     echo "Name of the sample ${split[0]}"
#     bash $DUP $PROJECTROOT ${split[0]}
#     echo " "
# done < <(tail -n +2 $EXPSUMMARYAL)
# 
# echo " "
# #Assess Fragment size
# 
# echo "Assess Fragment Size"
# while read line ; do
#      set $line
#      IFS=$','; split=($line); unset IFS;
#      if [[ $1 == "SampleName" ]]; then
#          continue
#      fi
#      echo "Name of the sample ${split[0]}"
#      bash $FRAG $PROJECTROOT ${split[0]}
#      echo " "
# done < <(tail -n +2 $EXPSUMMARYAL)
# 
# #Generate summary tables for both assessments
# echo "Summary Tables"
# Rscript $DUPSUM $PROJECTROOT
# Rscript $FRAGSUM $PROJECTROOT

#------------------4. Filter and convert---------------------------------------#
##minquality score to be filtered
# MINQUAL=2
# #Filter and convert & reproducibility assessment
# echo "Filter and convert"
# while read line ; do
#      set $line
#      IFS=$','; split=($line); unset IFS;
#      if [[ $1 == "SampleName" ]]; then
#          continue
#      fi
#      echo "Name of the sample ${split[0]}"
#      bash $FILTERCONV $MINQUAL $PROJECTROOT ${split[0]}
#      echo " "
# done < <(tail -n +2 $EXPSUMMARYAL)

#Create the correlation table to be plotted later 

#Rscript $REPREPRO $PROJECTROOT

#-------------------5. Normalization Igg and Spike in--------------------------#

##Spike-In alignment
cores=8
spikeinref=/project/genomes/Escherichia_coli/K12_DH10B/NCBI/2008-03-17/Sequence/Bowtie2Index
echo "Alignment to Spike in genome"
echo " "
mkdir -p $PROJECTROOT/alignment/sam/bowtie2_summary
mkdir -p $PROJECTROOT/alignment/bam
mkdir -p $PROJECTROOT/alignment/bed
mkdir -p $PROJECTROOT/alignment/bedgraph
while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     if [[ $1 == "SampleName" ]]; then
         continue
     fi
     echo "Name of the sample ${split[0]}"
     echo "Name of R1 ${split[1]}"
     echo "Name of R2 ${split[2]}"
     bash $SPIKEIN $RAWROOT/${split[1]} $RAWROOT/${split[2]} $PROJECTROOT \
     ${split[0]} $spikeinref
     echo " "
done < <(tail -n +2 $EXPSUMMARYAL)


#create a table summarizing the spike in alignments
echo "Spike-in summary"
Rscript $SPIKEINSUM $PROJECTROOT

##Igg substraction need to transform to bw files

echo "Convert to bigwig"

CHROMSIZES=/project/ChromGroup/Serkan_Project/cut_and_tag_analysis/tools/genome_annotations/hg38.chrom.sizes


while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     if [[ $1 == "SampleName" ]]; then
         continue
     fi
     bedGraphToBigWig $PROJECTROOT/alignment/bedgraph/${split[0]}_bowtie2.fragments.normalized.bedgraph \
     $CHROMSIZES $PROJECTROOT/alignment/bigwig/${split[0]}__bowtie2.fragments.normalized.bw
done < <(tail -n +2 $EXPSUMMARYAL)

echo "Igg substraction w/o spike in norm"




 
#------------------5. Peak Calling with SEACR----------------------------------#
