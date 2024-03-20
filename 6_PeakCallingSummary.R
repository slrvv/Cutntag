################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 6. Peak calling with SEACR                                                   #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: peak calling summary                                                #
################################################################################


#-------------------------Paths------------------------------------------------#
library(dplyr)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]

library(GenomicRanges)

#------------------------Sequencing depth--------------------------------------#
sampletable <- read.table(paste0(projPath, 
                                 "/experiment_summary_peaks.csv"),
                          header = T, sep = ",")

sampleList <- sampletable$SampleName

peakN = c()
peakWidth = c()
peakType = c("control", "top0.01")
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  if(histInfo[1] != "IgG"){
    for(type in peakType){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_seacr_", type, ".peaks.stringent.bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakType = type, Histone = histInfo[1], Replicate = histInfo[2]) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakType = type, Histone = histInfo[1], Replicate = histInfo[2])  %>% rbind(peakWidth, .)
    }
  }
}
peakN %>% select(Histone, Replicate, peakType, peakN)

histL = peakN$Histone
repL = c(1,2)
peakType = c("control", "top0.01")
peakOverlap = c()
for(type in peakType){
  for(hist in histL){
    overlap.gr = GRanges()
    for(rep in repL){
      file <- paste0(projPath, 
                     "/peakCalling/SEACR/",
                     hist,
                     "_",
                     rep,
                     "_seacr_",
                     type,
                     ".peaks.stringent.bed")

      if(!exists(file)){
        print(paste0(file, " doesn't exist"))
      }else{
        
        peakInfo = read.table(file, header = FALSE, fill = TRUE)
        peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
        if(length(overlap.gr) >0){
            overlap.gr = overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
        }else{
          overlap.gr = peakInfo.gr
          
      } 
      }
    }
    peakOverlap = data.frame(peakReprod = length(overlap.gr), Histone = hist, peakType = type) %>% rbind(peakOverlap, .)
    }
  }


peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "peakType")) %>% mutate(peakReprodRate = peakReprod/peakN * 100)
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)

write.table(peakReprod, paste0(projPath, 
                                 "alignment/summary/PeakCallingsummary.csv"),
                          header = T, row.names = F)
