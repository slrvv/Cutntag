################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling summary FRips Calculation                              #
################################################################################


#-------------------------Paths------------------------------------------------#
library(dplyr, quietly = TRUE)
library(chromVAR, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]
summaryPath <- args[2]


#------------------------Sequencing depth--------------------------------------#
sampletable <- read.table(summaryPath,
                          header = T, sep = ",")
alignResult <- read.table(paste0(projPath, 
                                 "/alignment/summary_seq_depth_all_experiments.txt"), 
                          header=T, sep = " ")
sampleList <- sampletable$SampleName
histL <- unique(sapply(strsplit(sampleList,"_"), `[`, 1))
repL = c(1,2)
bamDir = paste0(projPath, "/alignment/bam")
inPeakData = c()
## overlap with bam file to get count
for(hist in histL){
  for (rep in repL){
    peakRes = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_", rep, "_seacr_control.peaks.stringent.bed"), header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0(bamDir, "/", hist, "_", rep, ".mapped.bam")
    fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] %>% sum
    inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Histone = hist, Replicate = rep))
    
  }
    
}

frip = left_join(inPeakData, alignResult, by = c("Histone", "Replicate")) %>% mutate(frip = inPeakN/MappedFragNum_hg38 * 100)
frip %>% select(Histone, Replicate, SequencingDepth, MappedFragNum_hg38, AlignmentRate_hg38, FragInPeakNum = inPeakN, FRiPs = frip)

write.table(frip, paste0(projPath, 
                               "/alignment/summary_peak_calling_frips.txt"),
            row.names = F)