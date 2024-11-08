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
library(GenomicRanges, quiet = TRUE)
library(BRGenomics)
args <- commandArgs(trailingOnly=TRUE)
projPath <- args[1]
summaryPath <- args[2]

projPath <- "/project/ChromGroup_Seq_data/Celeste/2024_summer_HBD"
summaryPath <- "/project/ChromGroup_Seq_data/Celeste/2024_summer_HBD/experiment_summary_peaks_HBD.csv"

#------------------------Sequencing depth--------------------------------------#

sampletable <- read.table(summaryPath,
                          header = T, sep = ",")


sampleList <- sampletable$SampleName
sampleList
peakWidth = c()
peakN = c()
peakType = c("control", "top0.01", "control.rmDup", "top0.01.rmDup")
for(hist in sampleList){
  histInfo = unlist(strsplit(hist, "_", fixed = T))
  if(histInfo[1] != "IgG"){
    for(type in peakType){
      file <- paste0(projPath, "/peakCalling/SEACR/", hist, 
             "_seacr_", type, ".peaks.stringent.bed")
      print(file)
      if (file.info(file)$size != 0){
        peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, 
                                     "_seacr_", type, ".peaks.stringent.bed"),
                              header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
        print(nrow(peakInfo))
        peakN = data.frame(peakN = nrow(peakInfo), peakType = type, 
                           Histone = histInfo[1], 
                           Replicate = paste(histInfo[-1], collapse = "_")) %>% rbind(peakN, .)
        peakWidth = data.frame(width = peakInfo$width, 
                               peakType = type, 
                               Histone = histInfo[1], 
                               Replicate = paste(histInfo[-1], collapse = "_"))  %>% rbind(peakWidth, .)
      } else {
        peakN = data.frame(peakN = 0, peakType = type, 
                           Histone = histInfo[1], 
                           Replicate = paste(histInfo[-1], collapse = "_")) %>% rbind(peakN, .)
        peakWidth = data.frame(width = 0, 
                               peakType = type, 
                               Histone = histInfo[1], 
                               Replicate = paste(histInfo[-1], collapse = "_"))  %>% rbind(peakWidth, .)
      }
      
      
    }
  }
}


write.table(peakWidth, paste0(projPath, 
                               "/alignment/summary_peak_calling_width.txt"),
            row.names = F)

peakN %>% select(Histone, Replicate, peakType, peakN)

peakOverlap = c()

for(hist in sampleList){
  for(type in peakType){
    overlap.gr = GRanges()
    histInfo = unlist(strsplit(hist, "_", fixed = T))
    repL <- unique(peakN$Replicate[peakN$Histone == histInfo[1]])
    if(length(repL) > 1){
      for(rep in repL){
        file <- paste0(projPath,
                       "/peakCalling/SEACR/",
                       histInfo[1],
                       "_",
                       rep,
                       "_seacr_",
                       type,
                       ".peaks.stringent.bed")
        if (file.exists(file) & file.info(file)$size != 0){
          peakInfo = read.table(file, header = FALSE, fill = TRUE)
          peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
          if(length(overlap.gr) >0){
            overlap.gr = subsetByOverlaps(overlap.gr, peakInfo.gr)
            peakReprod = length(overlap.gr)
          }else{
            overlap.gr = peakInfo.gr
        
          }
          
        } else {
          peakReprod <- NaN
        }
        
      } 
    } else {
      peakReprod <- NaN
    }
    peakOverlap = data.frame(peakReprod = peakReprod, 
                             Histone = histInfo[1],
                             Replicate = histInfo[-1], 
                             peakType = type) %>% rbind(peakOverlap, .)
  }
}


peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "Replicate", "peakType")) %>% mutate(peakReprodRate = (peakReprod/peakN) * 100)
peakReprod
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)


write.table(peakReprod, paste0(projPath, 
                                 "/alignment/summary_peak_calling.txt"),
            row.names = F)
