################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 5. Normalization                                                             #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: substract igg                                                       #
################################################################################

#-------------------Paths------------------------------------------------------#
PROJECTROOT=$1
BWPATH=$PROJECTROOT/alignment/bigwig
FILE1=$2
IGG=$3
DUPREMOVE=$4
#----------------Igg substaction-----------------------------------------------#

#Spike in, no dedup
echo "Removing $3 from $2"
prun python bigwigCompare -b1 $BWPATH/${FILE1}_bowtie2.fragments.normalized.bw -b2 $BWPATH/${IGG}_bowtie2.fragments.normalized.bw  \
 --operation subtract  -p 8 -o $BWPATH/${FILE1}_bowtie2.fragments.normalized.substracted.igg.bw
 
#no spike in, no dedup
prun python bigwigCompare -b1 $BWPATH/${FILE1}_bowtie2.fragments.bw -b2 $BWPATH/${IGG}_bowtie2.fragments.bw \
 --operation subtract  -p 8 -o $BWPATH/${FILE1}_bowtie2.fragments.substracted.igg.bw

if [ "$DUPREMOVE" = true ]; then
  # spike in, dedup
  prun python bigwigCompare -b1 $BWPATH/${FILE1}_bowtie2.rmDup.fragments.normalized.bw -b2 $BWPATH/${IGG}_bowtie2.rmDup.fragments.normalized.bw \
 --operation subtract  -p 8 -o $BWPATH/${FILE1}_bowtie2.rmDup.fragments.normalized.substracted.igg.bw
 
  #no spike in, dedup
  prun python bigwigCompare -b1 $BWPATH/${FILE1}_bowtie2.rmDup.fragments.bw -b2 $BWPATH/${IGG}_bowtie2.rmDup.fragments.bw \
   --operation subtract  -p 8 -o $BWPATH/${FILE1}_bowtie2.rmDup.fragments.substracted.igg.bw
  
fi
