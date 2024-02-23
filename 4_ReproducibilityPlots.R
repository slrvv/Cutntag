################################################################################
#                                                                              #
# cut & tag pipeline R-loop project                                            #
#                                                                              #
# 4. Reproducibility assessment plot                                           #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
# purpose: create heatmap for reproducibility assessment plot                  #
################################################################################

library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)
library(corrplot)
library(reshape2)


#histones
repResult <- read.table("/project/ChromGroup/Serkan_Project/cut_and_tag_rloops/alignment/rep_reproducibility_all_experiments.txt",
                        header=T)

repResult_hist <- repResult[,c("H3K27ac_1", "H3K27ac_2", 
                               "H3K4me1_1", "H3K4me3_1", "H3K4me3_2",
                               "IgG_1", "IgG_2")]
repResult_hist <- repResult_hist[c("H3K27ac_1", "H3K27ac_2", 
                                   "H3K4me1_1", "H3K4me3_1", "H3K4me3_2",
                                   "IgG_1", "IgG_2"),]

repResult_hist <- as.matrix(repResult_hist)
corrplot(repResult_hist, method = "color",addCoef.col = "black", number.digits = 2, number.cex = 0.65)


repResult_prot <- repResult[,c("CR56_1", "CR56_2", "CR56_3",
                               "S173_1", "S173_2","Flag_1", "Flag_2",
                               "G4_1", "G4_2",
                               "IgG_1", "IgG_2")]

repResult_prot <- repResult_prot[c("CR56_1", "CR56_2", "CR56_3",
                                   "S173_1", "S173_2","Flag_1", "Flag_2",
                                   "G4_1", "G4_2",
                                   "IgG_1", "IgG_2"),]

repResult_prot <- as.matrix(repResult_prot)

corrplot(repResult_prot, method = "color",addCoef.col = "black", number.digits = 2, number.cex = 0.65)

repResult_KD <- repResult[,c( "Cont.siRNA.G4_1", "Cont.siRNA.G4_2",
                              "Cont.siRNA.Flag_1", "Cont.siRNA.Flag_2",
                             "SiAREG.G4_1", "SiAREG.G4_2",
                             "SiAREG.flag_1", "SiAREG.flag_2", 
                             "siMETTL3.G4_1", "siMETTL3.G4_2", 
                             "siMETTL3.flag_1", "siMETTL3.flag_2", 
                             "siTopo1.G4_1", "siTopo1.G4_2" , 
                             "siTopo1.flag_1", "siTopo1.flag_2", 
                             "IgG_1", "IgG_2")]

head(repResult_KD)
repResult_KD <- repResult_KD[c("Cont-siRNA-G4_1", "Cont-siRNA-G4_2",
                                "Cont-siRNA-Flag_1", "Cont-siRNA-Flag_2",
                                "SiAREG-G4_1", "SiAREG-G4_2",
                                "SiAREG-flag_1", "SiAREG-flag_2", 
                                "siMETTL3-G4_1", "siMETTL3-G4_2", 
                                "siMETTL3-flag_1", "siMETTL3-flag_2", 
                                "siTopo1-G4_1", "siTopo1-G4_2" , 
                                "siTopo1-flag_1", "siTopo1-flag_2", 
                                "IgG_1", "IgG_2"),]

repResult_KD <- as.matrix(repResult_KD)
corrplot(repResult_KD, method = "color",addCoef.col = "black", number.digits = 2, number.cex = 0.65)
