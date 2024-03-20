################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling                                                        #
################################################################################

#------------------------------Paths-------------------------------------------#

seacr="/project/ChromGroup/Serkan_Project/cut_and_tag_analysis/tools/SEACR/SEACR_1.3.sh"
projPath=$1
histControl=$2
histName=$3
mkdir -p $projPath/peakCalling/SEACR


#-------------------------Peak calling-----------------------------------------#

##Peak calling with SEACR ask if you should igg substract the bedgraph files.
echo "${histName}_bowtie2.fragments.normalized.bedgraph"
bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
$projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
non stringent $projPath/peakCalling/SEACR/${histName}_seacr_control.peaks

bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
0.01 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks

