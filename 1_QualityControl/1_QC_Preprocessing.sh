################################################################################
#                                                                              #
# cut & tag analysis pipeline R-loop identification project                    #
# 1. Preprocessing & QC                                                        #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: Preprocessing of the cut & tag data FASTQC and possible merging of  #
# replicates.                                                                  #
################################################################################


#----------------------------Tools path----------------------------------------#

#FastQC
FASTQC=/scratch/ngsvin/bin/FastQC/fastqc

#----------------------------Files path----------------------------------------#

#Paths are passed as command line arguments $1 $2 marks the order
#path where we save the results of FASTQC, this is a new folder. Dont use
#an already existing folder path.
OUTPUTPATH=$1
#path to the raw experiment data .fastq.gz
EXPPATH=$2

#--------------------------Preprocessing script--------------------------------#

#FASTQC for quality check

#make the path to save the FASTQC outputs
mkdir -p ${OUTPUTPATH}

$FASTQC -o $OUTPUTPATH -f fastq $EXPPATH

#Merging
#No merging for now, so we can assess the quality of replicates with IGV ask 
#Sarah what she wants. 
