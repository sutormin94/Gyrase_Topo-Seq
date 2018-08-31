###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Script visualizes GCSs association with different sets of transcription units (TUs):
#AG - all genes,
#GLE - genes with low transcription level,
#GHE - genes with high level of transcription,
#AO - all operons,
#OLE - operons characterized by low level of transcription,
#OHE - operons with high transcription level.
#Simultaneously the script is revealing association of GCSs with a range of TUs regions:
#Upstream - corresponds to USUS,
#TU start - corresponds to USGB,
#TU end - corresponds to GBDS,
#Downstream - corresponds to DSDS.
###############################################

#######
#Libraries to be imported (or installed).
#######

#install.packages("readxl")
#install.packages("ggplot2")
#install.packages("gridExtra")
#install.packages("forcats")
library(readxl)
library(ggplot2)
library(gridExtra)
library(forcats)

#######
#Variables to be defined.
#######

#Path to input excell file and sheet
excell_path="C:/Sutor/science/DNA-gyrase/Results/GCSs_sets_and_motifs/GCSs_association_with_US_GB_DS/GCSs_association_with_US_GB_DS.xlsx"
excell_sheet_TU_sets="For_R_DOOR_del_cor_all_data"
excell_sheet_rRNA="rRNA_data"
excell_sheet_rRNA_Cfx_RifCfx="Cfx_RifCfx_rRNA_data"
excell_sheet_score="Score_R"

Outpath_TU_sets_GCSs_num="C:/Sutor/science/DNA-gyrase/Results/GCSs_sets_and_motifs/GCSs_association_with_US_GB_DS/GCSs_association_with_TUs_USUS_USGB_GBDS_DSDS_test.png"
Outpath_rRNA_GCSs_num="C:/Sutor/science/DNA-gyrase/Results/GCSs_sets_and_motifs/GCSs_association_with_US_GB_DS/GCSs_association_with_rRNA_US_GB_DS_test.png"
Outpath_rRNA_GCSs_num_Cfx_RifCfx="C:/Sutor/science/DNA-gyrase/Results/GCSs_sets_and_motifs/GCSs_association_with_US_GB_DS/Cfx_RifCfx_GCSs_association_with_rRNA_US_GB_DS_test.png"
Outpath_TU_sets_GCSs_score="C:/Sutor/science/DNA-gyrase/Results/GCSs_sets_and_motifs/GCSs_association_with_US_GB_DS/GCSs_score_association_with_TUs_USUS_USGB_GBDS_DSDS_test.png"

#######
#Imports data for GCSs association with TUs sets.
#######

GCSs_assoc_num=read_excel(excell_path, sheet=excell_sheet_TU_sets)
GCSs_assoc_num <-data.frame(GCSs_assoc_num)

#######
#Plots.
#######

p <- ggplot(data = GCSs_assoc_num, aes(x=fct_inorder(GCSs_assoc_num$Set), y=GCSs_assoc_num$Cfx, fill=fct_inorder(GCSs_assoc_num$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Number of GCSs\n(normalized)") +
  scale_fill_brewer(name="Region", palette="GnBu") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p1 <- ggplot(data = GCSs_assoc_num, aes(x=fct_inorder(GCSs_assoc_num$Set), y=GCSs_assoc_num$RifCfx, fill=fct_inorder(GCSs_assoc_num$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Number of GCSs\n(normalized)") +
  scale_fill_brewer(name="Region", palette="YlGn") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p2 <- ggplot(data = GCSs_assoc_num, aes(x=fct_inorder(GCSs_assoc_num$Set), y=GCSs_assoc_num$Microcin, fill=fct_inorder(GCSs_assoc_num$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Number of GCSs\n(normalized)") +
  scale_fill_brewer(name="Region", palette="Reds") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p3 <- ggplot(data = GCSs_assoc_num, aes(x=fct_inorder(GCSs_assoc_num$Set), y=GCSs_assoc_num$Oxo, fill=fct_inorder(GCSs_assoc_num$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Number of GCSs\n(normalized)") +
  scale_fill_brewer(name="Region", palette="Blues") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p_sum <- grid.arrange(p, p1, p2, p3, nrow = 4)  
ggsave(Outpath_TU_sets_GCSs_num, units="in", p_sum, width=10, height=8, dpi=600)

#######
#Imports data for GCSs association with rRNA operons.
#######

Num_peaks_rRNA=read_excel(excell_path, sheet=excell_sheet_rRNA)
Num_peaks_rRNA <-data.frame(Num_peaks_rRNA)

#######
#Plots.
#######

p <- ggplot(data=Num_peaks_rRNA, aes(x=fct_inorder(Num_peaks_rRNA$Set) , y=Num_peaks_rRNA$Enrichment, fill=fct_inorder(Num_peaks_rRNA$Region))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Number of GCSs") +
  scale_fill_brewer(name="Region", palette="Set2") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p_sum_rRNA <- grid.arrange(p, nrow=1)
ggsave(Outpath_rRNA_GCSs_num, units="in", p_sum_rRNA, width=10, height=4, dpi=600)


#######
#Imports data for Cfx and RifCfx GCSs association with rRNA operons.
#######

Num_peaks_rRNA_CR=read_excel(excell_path, sheet=excell_sheet_rRNA_Cfx_RifCfx)
Num_peaks_rRNA_CR <-data.frame(Num_peaks_rRNA_CR)

#######
#Plots.
#######

p <- ggplot(data=Num_peaks_rRNA_CR, aes(x=fct_inorder(Num_peaks_rRNA_CR$Region) , y=Num_peaks_rRNA_CR$Enrichment, fill=fct_inorder(Num_peaks_rRNA_CR$Set))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=1, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Enrichment in the\nnumber of GCSs") +
  scale_fill_brewer(name="Condition:", palette="Set2") +
  theme(axis.line=element_line(size=1), axis.ticks=element_line(size=1), axis.text.x=element_text(size=20, color='black'),
        axis.text.y=element_text(size=23, color='black'), axis.title.y=element_text(size=25), axis.title.x=element_text(size=25),
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.2, "cm"),
        legend.text=element_text(size=25), legend.title=element_text(size=28,face = "bold"),
        legend.position=c(0.17,0.8), legend.key.size=unit(1, "cm"),
        plot.margin = margin(30, 0, 10, 0))

p_sum_rRNA <- grid.arrange(p, nrow=1)
ggsave(Outpath_rRNA_GCSs_num_Cfx_RifCfx, units="in", p_sum_rRNA, width=8, height=3.5, dpi=600)

#######
#Imports data for GCSs score association with TUs sets and compartments.
#######

Score=read_excel(excell_path, sheet=excell_sheet_score)
Score <-data.frame(Score)

#######
#Plots.
#######

p <- ggplot(data = Score, aes(x=fct_inorder(Score$Set), y=Score$Microcin, fill=fct_inorder(Score$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Average GCSs score") +
  scale_fill_brewer(name="Region", palette="Reds") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p1 <- ggplot(data = Score, aes(x=fct_inorder(Score$Set), y=Score$Cfx, fill=fct_inorder(Score$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Average GCSs score") +
  scale_fill_brewer(name="Region", palette="Greens") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p2 <- ggplot(data = Score, aes(x=fct_inorder(Score$Set), y=Score$Oxo, fill=fct_inorder(Score$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9)) +
  labs(x = NULL, y = "Average GCSs score") +
  scale_fill_brewer(name="Region", palette="Blues") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p_sum_score <- grid.arrange(p, p1, p2, nrow = 3) 
ggsave(Outpath_TU_sets_GCSs_score, units="in", p_sum_score, width=10, height=6.5, dpi=600)

