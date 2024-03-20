################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# processing of Laura's data                                                   #
# wetlab: Previous round (I don't know who it was)                             #
# purpose: Process the Flag and IgG data from the previous experimental round  #
################################################################################


#--------------------------Project Paths---------------------------------------#
#project folder
PROJECTROOT=/project/ChromGroup/Serkan_Project/cut_and_tag_rloops
#rootfolder to the raw experimental data
RAWROOT=/project/solexawork/pipelined/230302_SN435_7_AACJTW7M5/demux_20230303_SN435-AACJTW7M5_std-8d1/Project_OWL
PIPELINE=/project/ChromGroup/Serkan_Project/cut_and_tag_rloops/pipeline
#csv file with all of the experiments and a short id name

#csv file with all of the experiments but replicates are given in columns
#more suitable for the alignment step
EXPSUMMARY=$PROJECTROOT/previous_round_summary.csv



#-------------------------Module Paths-----------------------------------------#
#bash file that does the first step QC
QC=$PIPELINE/1_QC_Preprocessing.sh
#bash file for alignment
ALIGN=$PIPELINE/2_Alignment.sh
#bash file for spike in calibration
SPIKEIN=$PIPELINE/2bis_SpikeInAlignment_mm10.sh
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
#Substract Igg 
IGGSUB=$PIPELINE/5_IgGSubstract.sh
IGGSUBR=$PIPELINE/5_IgGSubstract.R
#Peak calling modules
PEAKS=$PIPELINE/6_PeakCalling.sh
PEAKSUMM=$PIPELINE/6_PeakCallingSummary.R
PEAKSFRIP=$PIPELINE/6_PeakCallingFrips.R


#---------------------------2. Alignment---------------------------------------#

# Alignment for each experiment, each experiment has 2 replicates EXPSUMMARYAL
# has one row per experiment and columns per replicate. split[0] is the name
# of the experiment split[1] is first replicate and split[2] is the second.
# The experiments were repeated by two different people which is marked by the _1
# or _2 in the name (split[0])
cores=8
ref=/project/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index_Large/genome
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
done < <(tail -n +2 $EXPSUMMARY)


#------------------4. Filter and convert---------------------------------------#
##minquality score to be filtered
MINQUAL=2
#Filter and convert & reproducibility assessment
echo "Filter and convert"
while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     if [[ $1 == "SampleName" ]]; then
         continue
     fi
     echo "Name of the sample ${split[0]}"
     bash $FILTERCONV $MINQUAL $PROJECTROOT ${split[0]}
     echo " "
done < <(tail -n +2 $EXPSUMMARY)

#-------------------5. Normalization Igg and Spike in--------------------------#

#Spike-In alignment
cores=8

echo "Alignment to Spike in genome"
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
     bash $SPIKEIN  $PROJECTROOT $RAWROOT/${split[1]} \
     $RAWROOT/${split[2]} ${split[0]}
     echo " "
done < <(tail -n +2 $EXPSUMMARY)

##Igg substraction need to transform to bw files

echo "Convert normalized to bigwig"

CHROMSIZES=/project/ChromGroup/Serkan_Project/cut_and_tag_analysis/tools/genome_annotations/hg38.chrom.sizes

while read line ; do
     set $line
     IFS=$','; split=($line); unset IFS;
     bedGraphToBigWig $PROJECTROOT/alignment/bedgraph/${split[0]}_bowtie2.fragments.normalized.bedgraph \
     $CHROMSIZES $PROJECTROOT/alignment/bigwig/${split[0]}_bowtie2.fragments.normalized.bw
done < <(tail -n +2 $EXPSUMMARY)

#convert Laura's Igg to bw


#bamCoverage -b /project/ChromGroup/Laura/bioinformatic_analyses/Cut_Tag_2303/2_bams/TR17-IgG-1_S1.bt2.filt.noDup.bam -o  $PROJECTROOT/alignment/bigwig/TR17-IgG-1_S1.bt2.filt.noDup.bw

echo "Igg substraction with spike in norm"


echo "Igg substraction with spike in norm"
while read line ; do
    set $line
    IFS=$','; split=($line); unset IFS;
    bash $IGGSUB $PROJECTROOT  ${split[0]}_bowtie2.fragments.normalized.bw  \
IgG_prev_bowtie2.fragments.normalized.bw
done < <(tail -n +2 $EXPSUMMARYPEAKS)

#some minus values appear in the tracks. set them to 0 using IGGSUBR
while read line ; do
  set $line
  IFS=$','; split=($line); unset IFS;
  Rscript $IGGSUBR $PROJECTROOT ${split[0]}_bowtie2.fragments.normalized.bw.substracted.igg.bw
done < <(tail -n +2 $EXPSUMMARY)


#------------------6. Peak Calling with SEACR----------------------------------#

# Call the peaks with SEACHR
echo "SEACR Peak Calling"

bash $PEAKS $PROJECTROOT  IgG_prev Flag_prev



