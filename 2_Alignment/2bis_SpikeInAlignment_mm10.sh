################################################################################
#                                                                              #
# cut & tag analysis pipeline R-loop identification project                    #
# 2bis. Spike-In Alignment                                                     #
#                                                                              #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: Alignment to Spike-In genome for Spike-In calibration (only needed) #
# if the experiment was performed with a Spike-In)                             #
################################################################################



#-----------------------------Paths--------------------------------------------#
#reference genome of the spike-in
spikeInRef=/project/genomes/mm10/Sequence/Bowtie2Index_Large/genome
#chromsize files for the reference genome of experiment (human)
chromSize=/project/ChromGroup/Serkan_Project/cut_and_tag_analysis/tools/genome_annotations/hg38.chrom.sizes
#root path of th whole project
projPath=$1
#fastqfiles for the replicates
rep1=$2
rep2=$3
#short name to id experiment
name=$4
#------------------------Spike-in calibration script---------------------------#
cores=8

echo "Align to mm10"

echo $projPath/alignment/sam/${name}_bowtie2_spikeIn.sam
## bowtie2-build path/to/Ecoli/fasta/Ecoli.fa /path/to/bowtie2Index/Ecoli
bowtie2  --local --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 $rep1 -2 $rep2 -S $projPath/alignment/sam/${name}_bowtie2_spikeIn.sam &>  $projPath/alignment/sam/bowtie2_summary/${name}_bowtie2_spikeIn.txt

echo "Calculate seq depth of spike in genome"
##calculate seq depth of spike in genome
seqDepthDouble=`samtools view -F 0x04 $projPath/alignment/sam/${name}_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${name}_bowtie2_spikeIn.seqDepth


echo "Rescale"
if [[ "$seqDepth" -gt "1" ]]; then

  mkdir -p $projPath/alignment/bedgraph
  
  scale_factor=`echo "10000 / $seqDepth" | bc -l`
  echo "Scaling factor for $name is: $scale_factor!"
  bedtools genomecov -bg -scale $scale_factor -i \
  $projPath/alignment/bed/${name}_bowtie2.fragments.bed -g $chromSize > \
  $projPath/alignment/bedgraph/${name}_bowtie2.fragments.normalized.bedgraph

fi

