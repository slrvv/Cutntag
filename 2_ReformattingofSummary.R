################################################################################
#                                                                              #
# cut&tag processing pipeline R loop project                                   # 
# purpose: transform experiment_summary.txt into a format more suitable        #
# for the aligner                                                              #
################################################################################

#-------------------------Paths and libraries----------------------------------#

metadatpath <- "/project/ChromGroup/Serkan_Project/cut_and_tag_rloops/experiment_summary.csv"

#-----------------------Script start-------------------------------------------#

metadata <- read.csv(metadatpath, sep = ",")
sample <- strsplit(metadata$Sample, "_")
sample
df <- data.frame(matrix(unlist(sample), nrow=length(sample), byrow=TRUE))

df$X4 <- apply( df[ , c("X1", "X2") ], 1, paste , collapse = "_" )
SampleName <- unique(df$X4)

FileNameRep1 <- c()
FileNameRep2 <- c()

for (element in SampleName){
  rep1<-paste(element, 1, sep="_")
  rep2<-paste(element, 2, sep="_")
  FileNameRep1 <- c(FileNameRep1, metadata[metadata$Sample == rep1, 2])
  FileNameRep2 <- c(FileNameRep2, metadata[metadata$Sample == rep2, 2])
  
}

metadata_formatted <- as.data.frame(cbind(SampleName, 
                                          FileNameRep1, 
                                          FileNameRep2))
write.csv(metadata_formatted, 
          "/project/ChromGroup/Serkan_Project/cut_and_tag_rloops/experiment_summary_align_formatted.csv",
          row.names = F,
          quote = F)
