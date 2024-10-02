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
dupRemove=$4
mkdir -p $projPath/peakCalling/SEACR


#-------------------------Peak calling-----------------------------------------#

echo "Peak calling on normalized experiments"

bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
$projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
non stringent $projPath/peakCalling/SEACR/${histName}_seacr_control.peaks

bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
0.01 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks

if [ "$dupRemove" = true ]; then
  echo "Peak calling on normalized and deduplicated experiments"
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph \
  $projPath/alignment/bedgraph/${histControl}_bowtie2.rmDup.fragments.normalized.bedgraph \
  non stringent $projPath/peakCalling/SEACR/${histName}_seacr_control.rmDup.peaks
  
  bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph \
  0.01 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.rmDup.peaks
fi
