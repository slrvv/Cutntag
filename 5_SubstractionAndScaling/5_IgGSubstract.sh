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
#----------------Igg substaction-----------------------------------------------#

prun python bigwigCompare -b1 $BWPATH/$FILE1_bowtie2.fragments.normalized.bw -b2 $BWPATH/$IGG \
 --operation subtract  -p 8 -o $BWPATH/${FILE1}_bowtie2.fragments.normalized.substracted.igg.bw