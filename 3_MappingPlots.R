################################################################################
#                                                                              #
# cut & tag analysis pipeline R-loop identification project                    #
# 3. Mapping Summary Plots                                                     #
# wetlab: Serkan Meydaner & Celeste Franconi                                   #
#                                                                              #
# purpose: R script to compute plots of summary and statistics of the mappings #
################################################################################

#-----------------------Paths and libraries------------------------------------#
library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

alignResultpath <- "/project/ChromGroup/Serkan_Project/cut_and_tag_rloops/alignment/summary_seq_depth_all_experiments.txt"

alignResult <- read.table(alignResultpath, header=T)

AlignmentRate_cont <- strsplit(alignResult$AlignmentRate_hg38,
                                         "%")
AlignmentRate_cont <- unlist(AlignmentRate_cont)
AlignmentRate_cont <- as.numeric(AlignmentRate_cont)

alignResult$AlignmentRate_cont <- AlignmentRate_cont


fig3A = alignResult %>% ggplot(aes(x = Histone,
                                   y = SequencingDepth/1000000, 
                                   fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete= T,begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Sequencing Depth per Million") +
  xlab("") + 
  ggtitle("A. Sequencing Depth")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.7, size=10),
        legend.title=element_text(size=14))

fig3A

fig3B = alignResult %>% ggplot(aes(x = Histone,
                                   y = MappedFragNum_hg38/1000000,
                                   fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Mapped Fragments per Million") +
  xlab("") +
  ggtitle("B. Alignable Fragment (hg38)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.7, size=10),
        legend.title=element_text(size=14))
fig3B

fig3C = alignResult %>% ggplot(aes(x = Histone, 
                                   y = AlignmentRate_cont, 
                                   fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = as.factor(Replicate)), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = T,begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (hg38)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.7, size=10),
        legend.title=element_text(size=14))
fig3C
ggarrange(fig3A, fig3B, fig3C, ncol = 2, nrow=2, 
          common.legend = TRUE, legend="bottom")
