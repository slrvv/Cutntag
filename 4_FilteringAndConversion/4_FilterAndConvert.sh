################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 4. Filter & Convert                                                          #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: module to filter the sam files and convert them to bam & bed        #
# for the analysis of the Cut&tag data. Paired-end sequencing and no Spike-In  #
################################################################################

#---------------------------Paths and tools------------------------------------#

minQualityScore=$1
projPath=$2
filename=$3
DUPREMOVE=$4
SAMPATH=$projPath/alignment/sam
SAMPATHDEDUP=$projPath/alignment/removeDuplicate/
BAMPATH=$projPath/alignment/bam
BEDPATH=$projPath/alignment/bed
BWPATH=$projPath/alignment/bigwig
BEDGRAPHPATH=$projPath/alignment/bedgraph
CHROMSIZES=/project/ChromGroup/Serkan_Project/cut_and_tag_analysis/tools/genome_annotations/hg38.chrom.sizes

#------------------Filter and convert the files--------------------------------#

echo "Filter and convert the sam files using samtools"

echo "Filter $filename"
samtools view -h -q $minQualityScore ${SAMPATH}/${filename}_bowtie2.sam > \
${SAMPATH}/${filename}_QS$minQualityScore.sam

## Filter and keep the mapped read pairs
echo "Create a bam file"
samtools view -bS -F 0x04 $SAMPATH/${filename}_QS${minQualityScore}.sam > \
$BAMPATH/$filename.mapped.bam

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
echo "Binned bed files for replicate reproducibility assessment"
binLen=500
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' \
$projPath/alignment/bed/${filename}_bowtie2.fragments.bed | \
sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' | \
sort -k1,1V -k2,2n  > \
$projPath/alignment/bed/${filename}_bowtie2.fragmentsCount.bin$binLen.bed

# Convert to bedgraph
bedtools genomecov -bg -i $BEDPATH/${filename}_bowtie2.fragments.bed \
				-g $CHROMSIZES > $BEDGRAPHPATH/${filename}_bowtie2.fragments.bedgraph

echo "Convert to bigwig"
mkdir -p $BWPATH

bedGraphToBigWig $BEDGRAPHPATH/${filename}_bowtie2.fragments.bedgraph \
$CHROMSIZES $BWPATH/${filename}_bowtie2.fragments.bw

#--------------Filter and conevert duplicate removed files---------------------#
#Convert and filter dup removed files

if [ "$DUPREMOVE" = true ]; then

  echo "Filter and convert the sam Deduplicated files using samtools"
  
  echo "Filter $filename"
  samtools view -h -q $minQualityScore ${SAMPATHDEDUP}/${filename}_bowtie2.sorted.rmDup.sam > \
  ${SAMPATHDEDUP}/${filename}_QS$minQualityScore.sorted.rmDup.sam
  
  ## Filter and keep the mapped read pairs
  echo "Create a bam file"
  samtools view -bS -F 0x04 $SAMPATHDEDUP/${filename}_QS${minQualityScore}.sorted.rmDup.sam > \
  $BAMPATH/$filename.mapped.rmDup.bam
  
  ## Convert into bed file format
  echo "Convert to bed file format"
  bedtools bamtobed -i $BAMPATH/${filename}.mapped.rmDup.bam -bedpe > \
  $BEDPATH/${filename}_bowtie2.rmDup.bed
  
  ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
  echo "Clean the bed file"
  awk '$1==$4 && $6-$2 < 1000 {print $0}' $BEDPATH/${filename}_bowtie2.rmDup.bed > \
  $BEDPATH/${filename}_bowtie2.rmDup.clean.bed
  
  ## Only extract the fragment related columns
  echo "Extract fragment related columns"
  cut -f 1,2,6 $BEDPATH/${filename}_bowtie2.rmDup.clean.bed | \
  sort -k1,1 -k2,2n -k3,3n  > $BEDPATH/${filename}_bowtie2.rmDup.fragments.bed
  
  ##bed file for replicate reproducibility assessment
  echo "Binned bed files for replicate reproducibility assessment"
  binLen=500
  awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' \
  $projPath/alignment/bed/${filename}_bowtie2.rmDup.fragments.bed | \
  sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' | \
  sort -k1,1V -k2,2n  > \
  $projPath/alignment/bed/${filename}_bowtie2.rmDup.fragmentsCount.bin$binLen.bed
  
  # Convert to bedgraph
  bedtools genomecov -bg -i $BEDPATH/${filename}_bowtie2.rmDup.fragments.bed \
  				-g $CHROMSIZES > $BEDGRAPHPATH/${filename}_bowtie2.rmDup.fragments.bedgraph
  
  echo "Convert to bigwig"
  mkdir -p $BWPATH
  
  bedGraphToBigWig $BEDGRAPHPATH/${filename}_bowtie2.rmDup.fragments.bedgraph \
  $CHROMSIZES $BWPATH/${filename}_bowtie2.rmDup.fragments.bw
fi
