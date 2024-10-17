#-------------------------------------------------------------------------------
# Author:  Sara Lopez Ruiz de Vargas
# Email:  lopez_s@molgen.mpg.de
#
# Date:    2024-04-12
#
# Script Name: cut_and_tag_pipeline_tutorial.sh
#
# Original Pipeline: https://yezhengstat.github.io/CUTTag_tutorial/#I_Introduction
#
# Script Description: This is a general use tutorial pipeline. Change the variables
# marked. Feel free to adapt the R plotting scripts so that your plots look their
# best.
#-------------------------------------------------------------------------------


#--------------------------Project Paths---------------------------------------#

#Reminder: In shell scripting variables are assigned using = without spaces
#Example: PROJECTROOT=/this/is/my/path will work.
#         PROJECTROOT = /this/is/my/path doesn't work.
# To access a variable we use '$',
# Example: echo $PROJECTROOT will print whatever is in that variable.
#          echo PROJECTROOT will throw an error or not work as expected.

#project folder
PROJECTROOT=/project/ChromGroup/Example 

#rootfolder to the raw experimental data (folder from SeqCore)

#Reminder: Raw data files are huge, there is no need to copy them into 
#your own directory, this leads to IT related problems.
RAWROOT=/project/solexawork/pipelined/230302_SN435_7_AACJTW7M5/demux_20230303_SN435-AACJTW7M5_std-8d1/Project_OWL/
#path to the folder with your pipeline
PIPELINE=$PROJECTROOT/pipeline
#csv file with all of the experiments and a short id name
EXPSUMMARY=$PROJECTROOT/experiment_summary_example.csv
#csv file with all of the experiments but replicates are given in columns
#more suitable for the alignment step
EXPSUMMARYAL=$PROJECTROOT/experiment_summary_align_formatted_example.csv
#csv file for peaks 
EXPSUMMARYPEAKS=$PROJECTROOT/experiment_summary_peaks_example.csv
DUPREMOVE=true #set to true if you want to have the duplicate removal files converted
IGGNEED=true #set to true if you want to remove IgG from bw files


#-------------------------Module Paths-----------------------------------------#
# No need to change anything from now on
#bash file that does the first step QC
QC=$PIPELINE/1_QualityControl/1_QC_Preprocessing.sh
#bash file for alignment
ALIGN=$PIPELINE/2_Alignment/2_Alignment.sh
#bash file for spike in calibration
SPIKEIN=$PIPELINE/2_Alignment/2bis_SpikeInAlignment.sh
#makes a summary of seq depth
SUMMARYAL=$PIPELINE/3_AssessAlignment/3_MappingSummary.R
#assess duplicates
DUP=$PIPELINE/3_AssessAlignment/3_DuplicateAssess.sh
#Fragment assessment
FRAG=$PIPELINE/3_AssessAlignment/3_FragAssess.sh
#summary of fragment and duplicates assessment
DUPSUM=$PIPELINE/3_AssessAlignment/3_DuplicateAssessSumary.R
FRAGSUM=$PIPELINE/3_AssessAlignment/3_FragAssessSumary.R
SPIKEINSUM=$PIPELINE/3_AssessAlignment/3_SpikeinSummary.R
#Filterandconvert
FILTERCONV=$PIPELINE/4_FilteringAndConversion/4_FilterAndConvert.sh
#Replicate reproducibility assessment
REPREPRO=$PIPELINE/4_FilteringAndConversion/4_ReplicateReproducibility.R
#Substract Igg 
IGGSUB=$PIPELINE/5_SubstractionAndScaling/5_IgGSubstract.sh
IGGSUBR=$PIPELINE/5_SubstractionAndScaling/5_IgGSubstract.R
#Peak calling modules
PEAKS=$PIPELINE/6_PeakCalling/6_PeakCalling.sh
PEAKSUMM=$PIPELINE/6_PeakCalling/6_PeakCallingSummary.R
PEAKSFRIP=$PIPELINE/6_PeakCalling/6_PeakCallingFrips.R

#-------------------------1. Quality Control-----------------------------------#

#Quality control each experiment. EXPSUMMARY is a table containing the file names
#for each experiment and a summary name. The FASTQC results are saved under
#the FastQC folder with subfolder name after the short experiment summary name
# split[0] is summary name split[1] is file name
echo "Quality Control with FASTQC"
echo " "
mkdir -p $PROJECTROOT/FastQCResults/
while read line ; do
    set $line
    IFS=$','; split=($line); unset IFS;
    echo "Name of experiment ${split[0]}"
    echo "Name of file ${split[1]}"
    bash ${QC} $PROJECTROOT/FastQCResults/${split[0]} $RAWROOT/${split[1]}
    echo " "
done < <(tail -n +2 $EXPSUMMARY)

#---------------------------2. Alignment---------------------------------------#

# Use bowtie2 : https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#end-to-end-alignment-example

#Alignment for each experiment, all experiments are paired so we have 2 files
#per experiment. EXPSUMMARYAL has one row per experiment and columns per file. split[0] is the name
#of the experiment split[1] is R1 and split[2] is R2.

cores=8 # Feel free to change this variable
# This is not only the path. You need to also state what name you use for
# your files. Example genome.1.bt2  genome.2.bt2 you should say "genome"
ref=/project/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index_Large/genome #you can also change this variable
echo "Alignment to hg38"
echo " "
mkdir -p $PROJECTROOT/alignment/sam/bowtie2_summary
mkdir -p $PROJECTROOT/alignment/bam
mkdir -p $PROJECTROOT/alignment/bed
mkdir -p $PROJECTROOT/alignment/bedgraph

while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     echo "Name of the sample ${split[0]}"
     echo "Name of R1 ${split[1]}"
     echo "Name of R2 ${split[2]}"
     bash $ALIGN $RAWROOT/${split[1]} $RAWROOT/${split[2]} $PROJECTROOT ${split[0]} $ref
     echo " "
done < <(tail -n +2 $EXPSUMMARYAL)

#-------------------3. Alignment Summary and assessment------------------------#

#Produce a table summarizing the sequencing depth of all of the alignments
echo "Make Alignment Summary"
Rscript $SUMMARYAL $PROJECTROOT $EXPSUMMARYAL
echo " "

#Create summary plots for the alignments using the 3_MappingPlots.R script

#Assess the number of duplicates

echo "Assess Duplicates"
while read line ; do
    set $line
    IFS=$','; split=($line); unset IFS;
    echo "Name of the sample ${split[0]}"
    bash $DUP $PROJECTROOT ${split[0]}
    echo " "
done < <(tail -n +2 $EXPSUMMARYAL)

echo " "
#Assess Fragment size

echo "Assess Fragment Size"
while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     echo "Name of the sample ${split[0]}"
     bash $FRAG $PROJECTROOT ${split[0]}
     echo " "
done < <(tail -n +2 $EXPSUMMARYAL)

#Generate summary tables for fragnment and duplicates statistics
echo "Summary Tables"
Rscript $DUPSUM $PROJECTROOT $EXPSUMMARYAL
Rscript $FRAGSUM $PROJECTROOT $EXPSUMMARYAL

#If you want to make plots to visualize fragment and duplicate statistics
#3_PlotsFragandDup.R

#------------------4. Filter and convert---------------------------------------#
#minquality score to be filtered
MINQUAL=2 ## Feel free to change this variable

#Filter and convert & reproducibility assessment

# Explanation of all of these formats: https://genome.ucsc.edu/FAQ/FAQformat.html
echo "Filter and convert"
while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     echo "Name of the sample ${split[0]}"
     bash $FILTERCONV $MINQUAL $PROJECTROOT ${split[0]} $DUPREMOVE
     echo " "
done < <(tail -n +2 $EXPSUMMARYAL)

#Create the correlation table between replicates
# We correlate the replicates for each experiment (remember that at this step no
# spike in or IgG substraction has been done) if your correlation with IgG is
# really high you might want to remove it.

#To study the reproducibility between replicates and across conditions,
#the genome is split into 500 bp bins, and a Pearson correlation of the 
#log2-transformed values of read counts in each bin is calculated between 
#replicate datasets. 
Rscript $REPREPRO $PROJECTROOT $EXPSUMMARYAL $DUPREMOVE

#Plot replicate matrices using 4_ReproducibilityPlots.R

#-------------------5. Normalization Igg and Spike in--------------------------#

#Spike-In alignment
cores=8 # feel free to change this variable
# Path to the Bowtie2 index that you are interested in.
# This is the one for e.coli
#if you have dna from another organism you can change refspike variable. Check
#which are available

#if you are using mouse:
refspike=/project/genomes/mm10/Sequence/Bowtie2Index/genome
#prebuild indexes https://benlangmead.github.io/aws-indexes/bowtie

#refspike=/project/genomes/Escherichia_coli/K12_DH10B/NCBI/2008-03-17/Sequence/Bowtie2Index/genome

#Human chromsizes
chromSize=$PROJECTROOT/pipeline/hg38.chrom.sizes


echo "Alignment to Spike in genome"
echo " "
while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     echo "Name of the sample ${split[0]}"
     echo "Name of R1 ${split[1]}"
     echo "Name of R2 ${split[2]}"
     bash $SPIKEIN  $PROJECTROOT $RAWROOT/${split[1]} \
     $RAWROOT/${split[2]} ${split[0]} $refspike $chromSize $DUPREMOVE
     echo " "
done < <(tail -n +2 $EXPSUMMARYAL)


#create a table summarizing the spike in alignments
echo "Spike-in summary"
Rscript $SPIKEINSUM $PROJECTROOT $EXPSUMMARYAL

#you can adapt 3_MappingPlots.R script to plot the summaries for spike-in

echo "Convert normalized to bigwig"

while read line ; do
 set $line
 IFS=$','; split=($line); unset IFS;
 bedGraphToBigWig $PROJECTROOT/alignment/bedgraph/${split[0]}_bowtie2.fragments.normalized.bedgraph \
 $chromSize $PROJECTROOT/alignment/bigwig/${split[0]}_bowtie2.fragments.normalized.bw
done < <(tail -n +2 $EXPSUMMARYAL)

# BW normalized, w/o Igg subs, with dedup
if [ "$DUPREMOVE" = true ]; then
  echo "Convert normalized and deduplicated to bigwig"
  while read line ; do
   set $line
   IFS=$','; split=($line); unset IFS;
   bedGraphToBigWig $PROJECTROOT/alignment/bedgraph/${split[0]}_bowtie2.rmDup.fragments.normalized.bedgraph \
   $chromSize $PROJECTROOT/alignment/bigwig/${split[0]}_bowtie2.rmDup.fragments.normalized.bw
  done < <(tail -n +2 $EXPSUMMARYAL)
fi

if [ "$IGGNEED" = true ]; then
  echo "Igg substraction"
  while read line ; do
       set $line
       IFS=$','; split=($line); unset IFS;
       bash $IGGSUB $PROJECTROOT ${split[0]}  \
       ${split[1]} $DUPREMOVE
  done < <(tail -n +2 $EXPSUMMARYPEAKS)

  # Set negative values to 0. You can choose to not do this part if you dont
  # have negative values.

  while read line ; do
       set $line
       IFS=$','; split=($line); unset IFS;
       Rscript $IGGSUBR $PROJECTROOT ${split[0]} $DUPREMOVE
  done < <(tail -n +2 $EXPSUMMARYPEAKS)
fi

#Check table with suffixes and their meaning on the README file

#------------------6. Peak Calling with SEACR----------------------------------#
# SEACR is run on the Spike in normalized bedgraphs only

# split[1] is Igg name and split[0] is sample name
# IgG_1 IgG_2
# If you have just one IgG replace split[1] with the name of your IgG file
# Ex: bash $PEAKS $PROJECTROOT IgG ${split[0]}
while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     echo "${split[0]}"
     bash $PEAKS $PROJECTROOT ${split[1]} ${split[0]} $DUPREMOVE
done < <(tail -n +2 $EXPSUMMARYPEAKS)

#Create summary statistics on the peaks called
Rscript $PEAKSUMM $PROJECTROOT $EXPSUMMARYPEAKS #Reproducibility and summary

# Frips: We calculate the fraction of reads in peaks (FRiPs) as a measure of signal-to-noise.
Rscript $PEAKSFRIP $PROJECTROOT $EXPSUMMARYPEAKS
# 6_PeakCallingSummaryPlot.R for visualization
