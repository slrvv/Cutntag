################################################################################
# cut & tag analysis pipeline R-loop identification project                    #
# 3. Aligment  assessment                                                      #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: Assessing Fragment size                                             #
################################################################################

#----------------------------Paths to project folders--------------------------#
### command line arguments
#projPath: project path folder containing all analysis on this particular
# experiment.
#name: short name abbreviation to identify experiment (i.e. H3K27ac)
projPath=$1
name=$2


#---------------------------Frag Size Script-----------------------------------#

#Technically one shouldn't remove duplicates in cutntag data. But we will assess
#whether in this case it makes sense or not.

mkdir -p $projPath/alignment/sam/fragmentLen

## Extract the 9th column from the alignment sam file which is the fragment length
samtools view -F 0x04 $projPath/alignment/sam/${name}_bowtie2.sam | \
 awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | \
 sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > \
 $projPath/alignment/sam/fragmentLen/${name}_fragmentLen.txt

