################################################################################
#                                                                              #
# cut & tag analysis pipeline R-loop identification project                    #
# 3. Frag Assess                                                               #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: R script to compute plots of the frag size and dups                 #
################################################################################

library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)

projectPath <- "/project/ChromGroup/Example" #fill this variable with the same path as PROJECTROOT in 
# the cut_and_tag_pipeline.sh script

dupResultPath <- paste0(projectPath, "/alignment/summary_duplicates_all_experiments.txt")

dupResult <- read.table(dupResultPath,
                        header=T)

DuplicationRate_cont <- strsplit(dupResult$DuplicationRate,
                               "%")
DuplicationRate_cont <- unlist(DuplicationRate_cont)

DuplicationRate_cont <- as.numeric(DuplicationRate_cont)
dupResult$Replicate <- as.factor(dupResult$Replicate)
dupResult$DuplicationRate_cont <- DuplicationRate_cont

fig4A = dupResult %>% ggplot(aes(x = Experiment, y = DuplicationRate_cont, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Duplication Rate (*100%)") +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.7, size=10),
        legend.title=element_text(size=14))

fig4A
fig4B = dupResult %>% ggplot(aes(x = Histone, y = EstimatedLibrarySize, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Estimated Library Size") +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.7, size=10),
        legend.title=element_text(size=14))

fig4B
fig4C = dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("# of Unique Fragments") +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.7, size=10),
        legend.title=element_text(size=14))
fig4C
ggarrange(fig4A, fig4B, fig4C, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")


## Fragment size plots

fragLenPath <- paste0(projectPath, "/alignment/summary_fraglen_all_experiments.txt")

fragLen <- read.table(fragLenPath,
                      header=T, stringsAsFactors = T)
fragLen$Replicate <- as.factor(fragLen$Replicate)


fig5A = fragLen %>% ggplot(aes(x = sampleInfo,
                               y = fragLen,
                               weight = Weight, 
                               fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ylab("Fragment Length") +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.7, size=10),
                legend.title=element_text(size=14))
fig5A

fig5B = fragLen %>% ggplot(aes(x = fragLen, 
                               y = fragCount, 
                               color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(linewidth = 0.5) + facet_wrap(Histone~., ncol = 4, nrow= 4,  scales="free") +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))
fig5B

ggarrange(fig5A, fig5B, ncol = 2)


