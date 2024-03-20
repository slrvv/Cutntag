################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling summary FRips Calculation                              #
################################################################################


#-------------------------Paths------------------------------------------------#
library(dplyr)
library(chromVAR)
library(GenomicRanges)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]

#------------------------Sequencing depth--------------------------------------#
sampletable <- read.table(paste0(projPath, 
                                 "/experiment_summary_peaks.csv"),
                          header = T, sep = ",")

sampleList <- sampletable$SampleName
histL <- unique(sapply(strsplit(sampleList,"_"), `[`, 1))
repL = paste0("rep", 1:2)
bamDir = paste0(projPath, "/alignment/bam")
inPeakData = c()
## overlap with bam file to get count
for(hist in sampleList){
    peakRes = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_seacr_control.peaks.stringent.bed"), header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0(bamDir, "/", hist, ".mapped.bam")
    fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] %>% sum
    inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Histone = hist))
}

frip = left_join(inPeakData, alignResult, by = c("Histone")) %>% mutate(frip = inPeakN/MappedFragNum_hg38 * 100)
frip %>% select(Histone, SequencingDepth, MappedFragNum_hg38, AlignmentRate_hg38, FragInPeakNum = inPeakN, FRiPs = frip)

write.table(frip, paste0(projPath, 
                               "alignment/summary/PeakCallingfripssummary.csv"),
            header = T, row.names = F)