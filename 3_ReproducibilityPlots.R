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


repResult <- read.table("/project/ChromGroup/Serkan_Project/cut_and_tag_rloops/alignment/rep_reproducibility_all_experiments.txt",
                        header=T)

repResult <- as.matrix(repResult)
meltedrepResult <- melt(repResult)
replicates <- c()

nameList <- sapply(strsplit(samplenames,"_"), `[`, 1)
frames <- cbind(nameList, samplenames)
frames
p2 = ggplot(data=repResult) +
  geom_raster(aes(x=Var1, y=Var2, fill=value)) +
  scale_fill_gradient2(low="blue", high="red", na.value="black", name="") +
  geom_rect(data=frames, size=1, fill=NA, colour="black",
            aes(xmin=Var1 - 0.5, xmax=Var1 + 0.5, ymin=Var2 - 0.5, ymax=Var2 + 0.5)) +
  labs(title="Example 2")

corrplot(repResult, method = "color",addCoef.col = "black", number.digits = 2, number.cex = 0.65)

corrplot(repResult, method = "color", outline = T, addgrid.col = "darkgray",
         order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,
         cl.pos = "b", tl.col = "indianred4", tl.cex = 1,
         cl.cex = 1, addCoef.col = "black", number.digits = 2,
         number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))
