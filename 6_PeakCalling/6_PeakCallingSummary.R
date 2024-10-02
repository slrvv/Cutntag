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
library(GenomicRanges)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]
summaryPath <- args[2]

#------------------------Sequencing depth--------------------------------------#
sampletable <- read.table(summaryPath,
                          header = T, sep = ",")

sampleList <- sampletable$SampleName

peakN = c()
peakWidth = c()
peakType = c("control", "top0.01", "control.rmDup", "top0.01.rmDup")
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

write.table(peakWidth, paste0(projPath, 
                               "/alignment/summary_peak_calling_width.txt"),
            row.names = F)

peakN %>% select(Histone, Replicate, peakType, peakN)

histL = unique(peakN$Histone)
repL = c(1,2)
peakType = c("control", "top0.01", "control.rmDup", "top0.01.rmDup")
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
        
        peakInfo = read.table(file, header = FALSE, fill = TRUE)
        peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
        if(length(overlap.gr) >0){
            overlap.gr = intersect(overlap.gr, peakInfo.gr)
        }else{
          overlap.gr = peakInfo.gr
          
        } 
    }
    peakOverlap = data.frame(peakReprod = length(overlap.gr), Histone = hist, peakType = type) %>% rbind(peakOverlap, .)
  }
}


peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "peakType")) %>% mutate(peakReprodRate = (peakReprod/peakN) * 100)
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)

write.table(peakReprod, paste0(projPath, 
                                 "/alignment/summary_peak_calling.txt"),
            row.names = F)
