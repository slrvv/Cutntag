################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# make igg merged                                                              #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: This script makes an IgG merged track for the correlation plot      #
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
EXPSUMMARYAL=$PROJECTROOT/experiment_summary_Latest.csv
EXPSUMMARYPEAKS=$PROJECTROOT/experiment_summary_peaks.csv


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
#Substract Igg 
IGGSUB=$PIPELINE/5_IgGSubstract.sh
IGGSUBR=$PIPELINE/5_IgGSubstract.R
#Peak calling modules
PEAKS=$PIPELINE/6_PeakCalling.sh
PEAKSUMM=$PIPELINE/6_PeakCallingSummary.R
PEAKSFRIP=$PIPELINE/6_PeakCallingFrips.R
SAMPATH=$PROJECTROOT/alignment/sam
BAMPATH=$PROJECTROOT/alignment/bam
BEDPATH=$PROJECTROOT/alignment/bed
BWPATH=$PROJECTROOT/alignment/bigwig
BEDGRAPHPATH=$PROJECTROOT/alignment/bedgraph
CHROMSIZES=/project/ChromGroup/Serkan_Project/cut_and_tag_analysis/tools/genome_annotations/hg38.chrom.sizes
filename=IgG_merged

#---------------------------Make iGg merged------------------------------------#


PROJECTROOT=/project/ChromGroup/Serkan_Project/cut_and_tag_rloops
#merge from bam
samtools merge $PROJECTROOT/alignment/bam/IgG_merged.mapped.bam \
$PROJECTROOT/alignment/bam/IgG_1.mapped.bam \
$PROJECTROOT/alignment/bam/IgG_2.mapped.bam

#file conversion and filtering

## Convert into bed file format
echo "Convert to bed file formant"
bedtools bamtobed -i $BAMPATH/${filename}.mapped.bam -bedpe > \
$BEDPATH/${filename}_bowtie2.bed

## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
echo "Clean the bed file"
awk '$1==$4 && $6-$2 < 1000 {print $0}' $BEDPATH/${filename}_bowtie2.bed > \
$BEDPATH/${filename}_bowtie2.clean.bed

## Only extract the fragment related columns
echo "Extract fragment related columns"
cut -f 1,2,6 $BEDPATH/${filename}_bowtie2.clean.bed | \
sort -k1,1 -k2,2n -k3,3n  > $BEDPATH/${filename}_bowtie2.fragments.bed

##bed file for replicate reproducibility assessment

# Convert to bedgraph
bedtools genomecov -bg -i $BEDPATH/${filename}_bowtie2.fragments.bed \
				-g $CHROMSIZES > $BEDGRAPHPATH/${filename}_bowtie2.fragments.bedgraph
				
#scale 



bedGraphToBigWig $BEDGRAPHPATH/${filename}_bowtie2.fragments.bedgraph \
$CHROMSIZES $BWPATH/${filename}_bowtie2.fragments.bw





