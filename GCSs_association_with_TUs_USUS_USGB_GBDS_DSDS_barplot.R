###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Script vizualises GCSs association with different sets of transcription units (TUs):
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
excell_path="C:/Sutor/science/DNA-gyrase/Results/Final_data_2/GCSS_association_with_US_GB_DS/DOOR_TUs/Deletions_corrected/Peaks_assoc_with_us_ds_gb.xlsx"
excell_sheet="For_R_DOOR_del_cor_all_data"

Outpath="C:/Sutor/science/DNA-gyrase/reports and burocraty/Papers/PICTURES/GCSs_association_with_US_DS_GB/test.png"

#######
#Imports data.
#######
GCSs_assoc_num=read_excel(excell_path, sheet=excell_sheet)

GCSs_assoc_num <-data.frame(GCSs_assoc_num)

#######
#Plots.
#######

p <- ggplot(data = GCSs_assoc_num, aes(x=fct_inorder(GCSs_assoc_num$Set), y=GCSs_assoc_num$Cfx, fill=fct_inorder(GCSs_assoc_num$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9))+
  labs(x = NULL, y = "Number of GCSs\n(normalized)") +
  scale_fill_brewer(name="Region", palette="GnBu") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p1 <- ggplot(data = GCSs_assoc_num, aes(x=fct_inorder(GCSs_assoc_num$Set), y=GCSs_assoc_num$RifCfx, fill=fct_inorder(GCSs_assoc_num$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9))+
  labs(x = NULL, y = "Number of GCSs\n(normalized)") +
  scale_fill_brewer(name="Region", palette="YlGn") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p2 <- ggplot(data = GCSs_assoc_num, aes(x=fct_inorder(GCSs_assoc_num$Set), y=GCSs_assoc_num$Microcin, fill=fct_inorder(GCSs_assoc_num$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9))+
  labs(x = NULL, y = "Number of GCSs\n(normalized)") +
  scale_fill_brewer(name="Region", palette="Reds") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p3 <- ggplot(data = GCSs_assoc_num, aes(x=fct_inorder(GCSs_assoc_num$Set), y=GCSs_assoc_num$Oxo, fill=fct_inorder(GCSs_assoc_num$Position))) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9))+
  labs(x = NULL, y = "Number of GCSs\n(normalized)") +
  scale_fill_brewer(name="Region", palette="Blues") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18, color="black"),
        axis.text.y=element_text(size=14), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"),
        legend.text=element_text(size=15), legend.title=element_text(size=17,face = "bold"))

p_sum <- grid.arrange(p, p1, p2, p3, nrow = 4)  

ggsave(Outpath, units="in", p_sum, width=10, height=8, dpi=600)

Num_peaks_rRNA=read_excel("C:/Sutor/science/DNA-gyrase/Results/Final_data_2/GCS_assocoation_with_highly_transcribed_operons/GCSs_association_with rRNA_operons_for_R.xlsx", sheet = "For_R")

Num_peaks_rRNA <-data.frame(Num_peaks_rRNA)
Num_peaks_rRNA$Set
Num_peaks_rRNA$Enrichment
Num_peaks_rRNA$Region

#rownames(Num_peaks) <- c("rRNA A", "rRNA B", "rRNA C", "rRNA D", "rRNA E", "rRNA G", "rRNA H", "Name", "Antibiotic")
#indata

p <- ggplot(data = Num_peaks_rRNA, aes(x = Num_peaks_rRNA$Set , y=Num_peaks_rRNA$Enrichment, fill= Num_peaks_rRNA$Region)) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9))+
  labs(x = NULL, y = "Number of GCSs") +
  scale_fill_brewer(name="Region", palette="Set2") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=11), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"))

p_sum <- grid.arrange(p, nrow = 1)  


Score=read_excel("C:/Sutor/science/DNA-gyrase/Results/Final_data_1/Association_with_US_IG_DS/Peaks_assoc_with_us_ds_gb.xlsx", sheet = "Score_R")

Score <-data.frame(Score)
Score$Set

#rownames(Num_peaks) <- c("rRNA A", "rRNA B", "rRNA C", "rRNA D", "rRNA E", "rRNA G", "rRNA H", "Name", "Antibiotic")
#indata

p <- ggplot(data = Score, aes(x = Score$Set , y=Score$Microcin, fill= Score$Position)) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9))+
  labs(x = NULL, y = "Average GCSs score") +
  scale_fill_brewer(name="Region", palette="Reds") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=11), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"))

p1 <- ggplot(data = Score, aes(x = Score$Set , y=Score$Cfx, fill= Score$Position)) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9))+
  labs(x = NULL, y = "Average GCSs score") +
  scale_fill_brewer(name="Region", palette="Greens") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=11), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"))
p2 <- ggplot(data = Score, aes(x = Score$Set , y=Score$Oxo, fill= Score$Position)) +
  geom_bar(stat = "identity", color="black", width=0.9, size=0.5, position = position_dodge(width = 0.9))+
  labs(x = NULL, y = "Average GCSs score") +
  scale_fill_brewer(name="Region", palette="Blues") +
  theme(axis.line=element_line(size=0.6), axis.ticks=element_line(size=0.4), axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=11), axis.title.y=element_text(size=15), 
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.15, "cm"))

p_sum <- grid.arrange(p, p1, p2, nrow = 3) 

