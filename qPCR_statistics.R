###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Script vizualises qPCR data.
###############################################

#######
#Libraries to be imported (or installed).
#######

#install.packages("ggplot2")
library(ggplot2)

#######
#Variables to be defined.
#######

#Output path for plot.
Outpath <- 'C:/Sutor/science/DNA-gyrase/reports and burocraty/Papers/PICTURES/qPCR/test.png'

#######
#qPCR data (units). 
#For raw data (Ct and transformations of it) look at Summary_table.xlsx, sheet DS1.
#######

#Cfx
CFX_Mu<-c(32.82, 19.93, 15.65)
CFX_2394<-c(2.77, 4.02, 3.81)
CFX_3602<-c(1, 1, 1.01)
CFX_3594<-c(14.97, 19.01, 26.39)
#RifCfx
RifCFX_Mu<-c(24.95, 23.57, 11.63)
RifCFX_2394<-c(1.78, 2.01, 1.61)
RifCFX_3602<-c(1, 1, 1.01)
RifCFX_3594<-c(3.86, 7.52, 4.58)

samples_list<-list(CFX_Mu, CFX_2394, CFX_3602, CFX_3594, RifCFX_Mu, RifCFX_2394, RifCFX_3602, RifCFX_3594)

#######
#Returns mean for the list, standard deviation and standard error
#######

Stat_function<-function(x) {
  x_mean<-mean(x)
  x_sd<-sd(x)
  x_se<-x_sd/sqrt(length(x))
  return(c(x_mean, x_sd, x_se, length(x)))
}

samples_data<-t(data.frame(lapply(samples_list, Stat_function)))

samples_data<-data.frame(samples_data, c("MuSGS", "ccmH", "rRNA A US", "rRNA A DS", 
                                         "MuSGS", "ccmH", "rRNA A US", "rRNA A DS"), 
                                       c("Cfx","Cfx","Cfx","Cfx", 
                                         "RifCfx", "RifCfx", "RifCfx", "RifCfx"))


rownames(samples_data)<-c("Cfx_Mu", "Cfx_2394", "Cfx_3602", "Cfx_3594", 
                          "RifCfx_Mu", "RifCfx_2394", "RifCfx_3602", "RifCfx_3594")
colnames(samples_data)<-c("mean", "sd", "se", "len", "name", "antib_index")

samples_data

#######
#Calculates boundaries for ~0.95 confidence interval (mean +/-2*SE)
#######

limits <- aes(ymax = samples_data$mean + 2*samples_data$se,
              ymin = samples_data$mean - 2*samples_data$se)

#######
#Calculates t-test between correspondence sets of values from Cfx and RifCfx groups. 
#Prepares data for plotting.
#######

p_value_labels <- c(replicate(4, "p-value"), '','','','')

p_value_y <- c(39, 11, 8, 33, 0, 0, 0, 0)
p_value_x <- c(2, 1, 4, 3, 1, 2, 3, 4)

p_value_labels_ <- c(round(t.test(CFX_Mu, RifCFX_Mu)$p.value, digits=3),
                     round(t.test(CFX_2394, RifCFX_2394)$p.value, digits=3),
                     round(t.test(CFX_3602, RifCFX_3602)$p.value, digits=3),
                     round(t.test(CFX_3594, RifCFX_3594)$p.value, digits=3), '','','','')
p_value_y_ <- c(35, 7, 4, 29, 0, 0, 0, 0)
p_value_x_ <- c(2, 1, 4, 3, 1, 2, 3, 4)

#######
#Plotting.
#######

p <- ggplot(data = samples_data, aes(x = name, y = mean, fill = antib_index)) +
  geom_bar(stat = "identity", color="black", width=0.9, size=1, position = position_dodge(width = 0.9)) +
  geom_errorbar(limits, position = position_dodge(width = 0.9), width = 0.3, size=0.8) +
  labs(x = "Site", y = "Fold enrichment, units") +
  scale_fill_brewer(palette="Greens", name="Condition:") +
  geom_text(aes(x=p_value_x, y=p_value_y, label=p_value_labels), 
                size=8) +
  geom_text(aes(x=p_value_x_, y=p_value_y_, label=p_value_labels_), 
            size=8) +
  theme(axis.line=element_line(size=1), axis.ticks=element_line(size=1), axis.text.x=element_text(size=23, color='black'),
        axis.text.y=element_text(size=25, color='black'), axis.title.y=element_text(size=30), axis.title.x=element_text(size=30),
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.2, "cm"), 
        legend.position=c(0.85,0.5), legend.text=element_text(size=30), 
        legend.title=element_text(size=30), legend.key.size=unit(1.2, "cm"))+
  scale_y_continuous(breaks=seq(0, 40, 10), limits=c(0,40))
p

ggsave(Outpath, units="in", width=8, height=4, dpi=600)






