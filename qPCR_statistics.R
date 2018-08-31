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

#Output path for plots.
Outpath <- 'C:/Sutor/science/DNA-gyrase/reports and burocraty/Papers/PICTURES/qPCR/qPCR_Cfx_RifCfx_test.png'
Outpath_qpcr_cfx <- 'C:/Sutor/science/DNA-gyrase/reports and burocraty/Papers/PICTURES/qPCR/qPCR_Cfx_test.png'
Outpath_cov_cfx <- 'C:/Sutor/science/DNA-gyrase/reports and burocraty/Papers/PICTURES/qPCR/cov_Cfx_test.png'

#######
#qPCR data (units). 
#For raw data (Ct and transformations of it) look at Summary_table.xlsx, sheet DS1.
#######

#Cfx Topo-qPCR fold enrichment
CFX_Mu<-c(32.82, 19.93, 15.65)
CFX_2394<-c(2.77, 4.02, 3.81)
CFX_3602<-c(1, 1, 1.01)
CFX_3594<-c(14.97, 19.01, 26.39)
#RifCfx Topo-qPCR fold enrichment
RifCFX_Mu<-c(24.95, 23.57, 11.63)
RifCFX_2394<-c(1.78, 2.01, 1.61)
RifCFX_3602<-c(1, 1, 1.01)
RifCFX_3594<-c(3.86, 7.52, 4.58)
#Cfx Topo-Seq mean coverage depth
Cov_CFX_Mu<-c(1203.0, 1783.508)
Cov_CFX_2394<-c(241.155, 189.264, 321.2)
Cov_CFX_3602<-c(195.245, 113.311, 210.311)
Cov_CFX_3594<-c(2271.133, 2205.693, 1466.08)

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

samples_data<-data.frame(samples_data, c("Mu SGS", "ccmH", "rRNA A US", "rRNA A DS", 
                                         "Mu SGS", "ccmH", "rRNA A US", "rRNA A DS"), 
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

p_value_label_y <- c(42, 13, 9, 36, 0, 0, 0, 0)
p_value_label_x <- c(2, 1, 4, 3, 1, 2, 3, 4)

p_value_values <- c(round(t.test(CFX_Mu, RifCFX_Mu)$p.value, digits=3),
                     round(t.test(CFX_2394, RifCFX_2394)$p.value, digits=3),
                     round(t.test(CFX_3602, RifCFX_3602)$p.value, digits=3),
                     round(t.test(CFX_3594, RifCFX_3594)$p.value, digits=3), '','','','')
p_value_value_y <- c(38, 9, 5, 32, 0, 0, 0, 0)
p_value_value_x <- c(2, 1, 4, 3, 1, 2, 3, 4)

#######
#Plotting barplot compares Cfx and RifCfx qPCR data.
#######

p <- ggplot(data = samples_data, aes(x = name, y = mean, fill = antib_index)) +
  geom_bar(stat = "identity", color="black", width=0.9, size=1, position = position_dodge(width = 0.9)) +
  geom_errorbar(limits, position = position_dodge(width = 0.9), width = 0.3, size=0.8) +
  labs(x = NULL, y = "Fold enrichment, \nunits") +
  scale_fill_brewer(palette="Set2", name="Condition:") +
  geom_text(aes(x=p_value_label_x, y=p_value_label_y, label=p_value_labels), 
                size=8) +
  geom_text(aes(x=p_value_value_x, y=p_value_value_y, label=p_value_values), 
            size=8) +
  theme(axis.line=element_line(size=1), axis.ticks=element_line(size=1), axis.text.x=element_text(size=20, color='black'),
        axis.text.y=element_text(size=25, color='black'), axis.title.y=element_text(size=25), axis.title.x=element_text(size=25),
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.2, "cm"), 
        legend.position=c(0.85,0.55), legend.text=element_text(size=25), 
        legend.title=element_text(size=28, face = "bold"), legend.key.size=unit(1, "cm"))+
  scale_y_continuous(breaks=seq(0, 45, 10), limits=c(0,42))
p

ggsave(Outpath, units="in", width=8, height=4, dpi=600)


#######
#Coverage data processing.
#######

#Cfx Topo-Seq mean coverage depth
Cov_CFX_Mu<-c(1203.0, 1783.508)
Cov_CFX_2394<-c(241.155, 189.264, 321.2)
Cov_CFX_3602<-c(195.245, 113.311, 210.311)
Cov_CFX_3594<-c(2271.133, 2205.693, 1466.08)

samples_list_cov<-list(Cov_CFX_Mu, Cov_CFX_2394, Cov_CFX_3602, Cov_CFX_3594)

#######
#Returns mean for the list, standard deviation and standard error
#######

Stat_function<-function(x) {
  x_mean<-mean(x)
  x_sd<-sd(x)
  x_se<-x_sd/sqrt(length(x))
  return(c(x_mean, x_sd, x_se, length(x)))
}

samples_data_cov<-t(data.frame(lapply(samples_list_cov, Stat_function)))

samples_data_cov<-data.frame(samples_data_cov, c("MuSGS", "ccmH", "rRNA A US", "rRNA A DS"), 
                         c("Cfx","Cfx","Cfx","Cfx"))


rownames(samples_data_cov)<-c("cov_Cfx_Mu", "cov_Cfx_2394", "cov_Cfx_3602", "cov_Cfx_3594")
colnames(samples_data_cov)<-c("mean", "sd", "se", "len", "name", "antib_index")

samples_data_cov

#######
#Calculates boundaries for ~0.95 confidence interval (mean +/-2*SE)
#######

limits_cov <- aes(ymax = samples_data_cov$mean + 2*samples_data_cov$se,
              ymin = samples_data_cov$mean - 2*samples_data_cov$se)

#######
#Plotting.
#######

p <- ggplot(data = samples_data_cov, aes(x = name, y = mean, fill = name)) +
  geom_bar(stat = "identity", color="black", width=0.95, size=1, position = position_dodge(width = 0.9)) +
  geom_errorbar(limits_cov, position = position_dodge(width = 0.95), width = 0.3, size=0.8) +
  labs(x = "Site", y = "Coverage depth, units") +
  scale_fill_brewer(palette="Blues", name="Condition:") +
  theme(axis.line=element_line(size=1), axis.ticks=element_line(size=1), axis.text.x=element_text(size=20, color='black'),
        axis.text.y=element_text(size=25, color='black'), axis.title.y=element_text(size=30), axis.title.x=element_text(size=30),
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.2, "cm"), 
        legend.position=c(0.2,0.8), legend.text=element_text(size=23), 
        legend.title=element_text(size=25), legend.key.size=unit(0.8, "cm"))+
  scale_y_continuous(breaks=seq(0, 2500, 400), limits=c(0,2500))
p

ggsave(Outpath_cov_cfx, units="in", width=8, height=8, dpi=600)


#######
#Cfx qPCR data processing.
#######

#Cfx Topo-qPCR fold enrichment
CFX_Mu<-c(32.82, 19.93, 15.65)
CFX_2394<-c(2.77, 4.02, 3.81)
CFX_3602<-c(1, 1, 1.01)
CFX_3594<-c(14.97, 19.01, 26.39)

samples_list_cov<-list(CFX_Mu, CFX_2394, CFX_3602, CFX_3594)

#######
#Returns mean for the list, standard deviation and standard error
#######

Stat_function<-function(x) {
  x_mean<-mean(x)
  x_sd<-sd(x)
  x_se<-x_sd/sqrt(length(x))
  return(c(x_mean, x_sd, x_se, length(x)))
}

samples_data_cov<-t(data.frame(lapply(samples_list_cov, Stat_function)))

samples_data_cov<-data.frame(samples_data_cov, c("MuSGS", "ccmH", "rRNA A US", "rRNA A DS"), 
                             c("Cfx","Cfx","Cfx","Cfx"))


rownames(samples_data_cov)<-c("cov_Cfx_Mu", "cov_Cfx_2394", "cov_Cfx_3602", "cov_Cfx_3594")
colnames(samples_data_cov)<-c("mean", "sd", "se", "len", "name", "antib_index")

samples_data_cov

#######
#Calculates boundaries for ~0.95 confidence interval (mean +/-2*SE)
#######

limits_cov <- aes(ymax = samples_data_cov$mean + 2*samples_data_cov$se,
                  ymin = samples_data_cov$mean - 2*samples_data_cov$se)

#######
#Plotting.
#######

p <- ggplot(data = samples_data_cov, aes(x = name, y = mean, fill = name)) +
  geom_bar(stat = "identity", color="black", width=0.95, size=1, position = position_dodge(width = 0.9)) +
  geom_errorbar(limits_cov, position = position_dodge(width = 0.9), width = 0.3, size=0.8) +
  labs(x = "Site", y = "Fold enrichment, units") +
  scale_fill_brewer(palette="Greens", name="Condition:") +
  theme(axis.line=element_line(size=1), axis.ticks=element_line(size=1), axis.text.x=element_text(size=20, color='black'),
        axis.text.y=element_text(size=25, color='black'), axis.title.y=element_text(size=30), axis.title.x=element_text(size=30),
        panel.background=element_rect(fill="white"), axis.ticks.length=unit(0.2, "cm"), 
        legend.position=c(0.85,0.8), legend.text=element_text(size=23), 
        legend.title=element_text(size=25), legend.key.size=unit(0.8, "cm"))+
  scale_y_continuous(breaks=seq(0, 40, 10), limits=c(0,35))
p

ggsave(Outpath_qpcr_cfx, units="in", width=8, height=8, dpi=600)




