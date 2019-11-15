###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes sets of trusted GCSs and analysis the distribution
#of GCSs trough the genome. Also it plots the distribution of other 
#values such as score, GC% and transcription.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.stats import binom
from scipy.stats import pearsonr

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Input data - GCSs, TAB.
path_to_GCSs_files={'Cfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\Cfx_10mkM_trusted_GCSs_h_s.txt",
                    'RifCfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\RifCfx_trusted_GCSs_h_s.txt",
                    'Micro': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\Micro_trusted_GCSs_h_s.txt",
                    'Oxo': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\Oxo_trusted_GCSs_h_s.txt"}                    

#Input data - score, WIG.
Score_data_path="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\E_coli_w3110_G_Mu_score.wig"

#Input data - GC, WIG.
GC_data_path="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\E_coli_w3110_Mu_GC_133bp.wig"

#Input data - transcription, TAB.
#Transcription_data_path="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_genes_expression.txt"
Transcription_data_path="F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\Expression_data\Exponential_phase_expression_TU_door_like_no_rRNA.txt"

#Input data - NAPs data in BroadPeak format.
NAPs_inpath={'IHFAB': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\IHFAB_sites_mid_exp_Prieto_W3110_Mu_SGS.BroadPeak",
             'H-NS': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\H-NS_sites_ME_Kahramanoglou_W3110_Mu_SGS.BroadPeak",
             'Fis': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Fis_sites_ME_Kahramanoglou_W3110_Mu_SGS.BroadPeak",
             'MatP': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MatP_sites_37deg_Nolivos_W3110_Mu_SGS.BroadPeak",
             'TopoIV' : "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\TopoIV_trusted_sites_Sayyed_W3110_Mu_SGS.BroadPeak",
             'MukB' : "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MukB_sites_37deg_Nolivos_W3110_Mu_SGS.BroadPeak"}

#Input data - chromosomal macrodomains.
Macrodomains_path="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Macrodomains_E_coli_w3110_edited.broadPeak"

#Input data - length of the genome.
Genome_length=4647454

#Output data: plot.
Plot_path_out="F:\Topo_data_new_expression\Exponential_phase_TUs_no_rRNA\\"
if not os.path.exists(Plot_path_out):
    os.makedirs(Plot_path_out)

#######
#Trusted GCSs data parsing.
#######

def trusted_GCSs_parsing(input_dict):
    GCSs_sets_dict={}
    for k, v in input_dict.items():
        ar=[]
        filein=open(v, 'r')
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in ['GCSs_coordinate']:
                ar.append(int(line[0]))
            else:
                continue
        GCSs_sets_dict[k]=ar
        print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(ar)))
    return GCSs_sets_dict

#######
#Parsing scores and/or GC - WIG file (optional step).
#######

def score_data_parser(inpath, param_name):
    param_file=open(inpath, 'r')
    ar=[]
    for line in param_file:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            ar.append(float(line[0]))
    param_file.close()
    print('Whole genome average ' + str(param_name) + ' : ' + str(sum(ar)/len(ar)))
    return ar     

#######
#Parsing transcription - TAB file (optional step).
#######

def transcription_data_parser(transcription_path):
    transcription_file=open(transcription_path, 'r')
    transcription=[]
    for i in range(4647454):
        transcription.append(0)
    for line in transcription_file:
        line=line.rstrip().split('\t')
        if line[0] not in ['GeneID', 'OperonID']:
            for j in range(int(line[3])-int(line[2])):
                transcription[int(line[2])+j]=float(line[5].replace(',','.'))
    transcription_file.close()
    print('Whole genome average transcription: ' + str(sum(transcription)/len(transcription)))
    return transcription


#######
#Parsing BroadPeak files (NAPs binding sites).
#Return array ready to be binned.
#######

def BroadPeak_to_bins(input_dict, bins):
    arrays_dict={}
    arrays_binned_dict={}
    for sample_name, path in input_dict.items():
        filein=open(path, 'r')
        ar=[]
        for i in range(4647454):
            ar.append(0)
        for line in filein:
            line=line.rstrip().split('\t')
            for j in range(int(line[2])-int(line[1])):
                ar[int(line[1])+j]=float(line[6])
        filein.close()
        print('Whole genome average ' + sample_name + ' : ' + str(sum(ar)/len(ar)))
        
        arrays_dict[sample_name]=ar
        ar_binned=bar_convert(ar, bins)
        arrays_binned_dict[sample_name]=ar_binned
    return arrays_dict, arrays_binned_dict


#######
#BroadPeak-formalized intervals parsing (macrodomains, etc) and filtering intervals that are not deleted.
#Return intervals.
#######

def broadpeak_pars(broadpeak_path, genome_len):
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    filein=open(broadpeak_path, 'r')
    ar_bins=[]   
    ar_names=[]
    ar_lens=[]
    for line in filein:
        line=line.rstrip().split('\t')
        int_start=int(line[1])
        int_end=int(line[2])
        MD_len=int_end-int_start
        #for j in range(len(deletions)):
        #    if deletions[j][1]>int_start>deletions[j][0]:
        #        MD_len+=int_start-deletions[j][1]
        #    elif deletions[j][1]>int_end>deletions[j][0]:
        #        MD_len+=deletions[j][0]-int_end
        #    elif int_end>deletions[j][1]>deletions[j][0]>int_start:
        #        MD_len+=int_start-int_end
        ar_bins.append(int_start) #Macrodomain starts which mark MDs bounadaries
        ar_names.append(line[3]) #Macrodomains names
        ar_lens.append(MD_len) #Macrodomains length corrected for deletions
    ar_bins.append(int_end)
    print("Number of chromosomal macrodomains: " + str(len(ar_names)))
    filein.close()
    Macrodomains_descr={"Bins": ar_bins, "Names": ar_names, "Lengths": ar_lens}
    return Macrodomains_descr

#######
#Convert score/transcription/GC data to histogram.
#######

def bar_convert(ar, bins):
    bar_ar=[]
    for i in range(len(bins)-1):
        if bins[i]>=0:
            bar_ar.append(sum(ar[int(bins[i]):int(bins[i+1])])/(bins[i+1]-bins[i]))
        else:
            bar_ar.append(sum(ar[int(bins[i]):-1]+ar[0:int(bins[i+1])])/(bins[i+1]-bins[i]))
    return bar_ar

#######
#Wrapps all the nonGCSs data together and prepare it for plotting.
#######

def Prepare_non_GCSs_data(score_path, GC_path, transcription_path, bins):
    Data_dict={}
    Score_data=score_data_parser(score_path, 'score')
    GC_data=score_data_parser(GC_path, 'GC')
    Trancription_data=transcription_data_parser(transcription_path)
    Data_dict['Score']=bar_convert(Score_data, bins)
    Data_dict['GC']=bar_convert(GC_data, bins)
    Data_dict['Transcription']=bar_convert(Trancription_data, bins)
    return Data_dict


#######
#Plotting: distribution throughout the genome, bins are evenly distributed.
#######

def Plot_the_distribution(GSCs_data, Non_GCSs_data, bins, genome_len, path_out):
    #Parameters
    ticks1=[0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000]
    xticknames1=['', '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500']
    colors=['#7FCE79', '#BAE85C', '#ff878b', '#8991ff', '#ac5eff', '#50b3ff', '#ffd75e']
    plot_names=['plot1', 'plot2', 'plot3', 'plot4', 'plot5', 'plot6', 'plot7']
    Y_labels=['Cfx GCSs', 'RifCfx GCSs', 'Micro GCSs', 'Oxo GCSs', 'Score', 'GC%', 'Transcription\nlevel']
    yter=1592477
    yori=3711828
    #GCSs data plotting.
    fig, plot_names=plt.subplots(7,1,figsize=(11,15), dpi=100)
    i=0
    Histo_comp_dict={} #Will contain computed histogramm data (bins and values)
    for key, value in GSCs_data.items():
        plot_names[i].set_xlim(0, genome_len)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticks([yter, yori], minor=True)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        conf_interval=[binom.interval(0.999, len(value), 1/10)[0], binom.interval(0.999, len(value), 1/10)[1]]
        plot_names[i].set_yticks(conf_interval, minor=True)
        plot_names[i].yaxis.grid(True, which='minor', linewidth=0.4, linestyle='--', color='black')
        plot_names[i].fill_between(bins, conf_interval[0], conf_interval[1], facecolor='grey', alpha=0.3)
        plot_names[i].locator_params(axis='y', nbins=6)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        Histo_comp_dict[key]=plot_names[i].hist(value, bins, facecolor=colors[i], alpha=0.7, linewidth=1, edgecolor='black') #Plot histo and save computed histogramm data (bins and values)
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(Y_labels[i], size=22, labelpad=8, rotation=90)
        i+=1
    #Score, GC, Transcription plotting.
    bin_width=int(bins[1])
    position=bins[:-1]
    print(len(position))
    for key, value in Non_GCSs_data.items():
        len(value)
        plot_names[i].set_xlim(0, genome_len)
        if key=="GC":
            plot_names[i].set_ylim(45, max(value)+2)
        elif key=="Score":
            plot_names[i].set_ylim(min(value)-0.2, -1.5)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticks([yter, yori], minor=True)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        plot_names[i].locator_params(axis='y', nbins=6)
        plot_names[i].bar(position, value, bin_width, color=colors[i], linewidth=1, edgecolor='black', align='edge')
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(Y_labels[i], size=22, labelpad=8, rotation=90) 
        i+=1
    plt.tight_layout()
    fig.savefig(path_out+"GCSs_num_score_GC133_transcription_distrib_thr_genome.png", figsize=(11,15), dpi=400) 
    
    #GCSs data plotting for Cfx, Micro, and Oxo only.
    GSCs_data_main={'Cfx': GSCs_data['Cfx'], 'Micro': GSCs_data['Micro'], 'Oxo': GSCs_data['Oxo']}
    Y_labels_main=['Cfx GCSs', 'Micro GCSs', 'Oxo GCSs']
    colors_main=['#7FCE79', '#ff878b', '#8991ff']
    fig, plot_names=plt.subplots(3,1,figsize=(11,7), dpi=100)
    i=0
    for key, value in GSCs_data_main.items():
        plot_names[i].set_xlim(0, genome_len)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticks([yter, yori], minor=True)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        conf_interval=[binom.interval(0.999, len(value), 1/10)[0], binom.interval(0.999, len(value), 1/10)[1]]
        plot_names[i].set_yticks(conf_interval, minor=True)
        plot_names[i].yaxis.grid(True, which='minor', linewidth=0.4, linestyle='--', color='black')
        plot_names[i].fill_between(bins, conf_interval[0], conf_interval[1], facecolor='grey', alpha=0.3)
        plot_names[i].locator_params(axis='y', nbins=6)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        plot_names[i].hist(value, bins, facecolor=colors_main[i], alpha=0.7, linewidth=1, edgecolor='black') #Plot histo and save computed histogramm data (bins and values)
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(Y_labels_main[i], size=22, labelpad=8, rotation=90)
        i+=1    
    plt.tight_layout()
    fig.savefig(path_out+"GCSs_number_Cfx_Micro_Oxo_distrib_thr_genome.png", figsize=(11,7), dpi=400) 
    return Histo_comp_dict


#######
#Plotting: NAPs distribution throughout the genome, bins are evenly distributed.
#######

def Plot_NAPs_distribution(Some_data, bins, path_out, genome_len, binning_type):
    #Parameters
    ticks1=[0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000]
    xticknames1=['', '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500']
    colors=['#7FCE79', '#BAE85C', '#ff878b', '#8991ff', '#ac5eff', '#50b3ff']
    plot_names=['plot1', 'plot2', 'plot3', 'plot4', 'plot5', 'plot6']
    yter=1592477
    yori=3711828
    
    #Prepare binning for plotting.
    if binning_type=='even':
        bin_width=int(bins[1])
        position=bins[:-1]
        print("Number of bins: " + str(len(position)))
    
    elif binning_type=='MDs':
        print("MDs binning: " + str(bins))
        position=[0]+bins[1:]
        print("Bars positions: " + str(position))
        bin_width=[]
        position_bw=position+[genome_len]
        for j in range(len(position_bw)-1):
            bin_width.append(position_bw[j+1]-position_bw[j])
        print("Bins width: " + str(bin_width))
        position_centre=[]
        for j in range(len(bin_width)):
            position_centre.append(int(position[j]+(bin_width[j]/2)))
        position_centre=[0]+position_centre+[genome_len]      
    
    #Data plotting.
    fig, plot_names=plt.subplots(6,1,figsize=(11,13), dpi=100)
    i=0

    for key, value in Some_data.items():
        if binning_type=='even':
            print("Number of " + key + " data bins: " + str(len(value)))
        elif binning_type=='MDs':
            value.append(value[0])
            print("Number of " + key + " data bins: " + str(len(value)))
        plot_names[i].set_xlim(0, genome_len)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticks([yter, yori], minor=True)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        plot_names[i].locator_params(axis='y', nbins=6)
        plot_names[i].bar(position, value, bin_width, color=colors[i], linewidth=1, edgecolor='black', align='edge')
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(key, size=22, labelpad=8, rotation=90) 
        i+=1
    plt.tight_layout()
    fig.savefig(path_out, figsize=(11,13), dpi=400)  
    plt.close()
    return 

#######
#Plotting: distribution throughout the genome, macrodomains are as bins.
#######

def Plot_the_distribution_MDs(GSCs_data, Non_GCSs_data, bins, MDs_lengths, genome_len, path_out):
    #Parameters
    ticks1=[0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000]
    xticknames1=['', '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500']
    colors=['#7FCE79', '#BAE85C', '#ff878b', '#8991ff', '#ac5eff', '#50b3ff', '#ffd75e']
    plot_names=['plot1', 'plot2', 'plot3', 'plot4', 'plot5', 'plot6', 'plot7']
    Y_labels=['Cfx GCSs', 'RifCfx GCSs', 'Micro GCSs', 'Oxo GCSs', 'Score', 'GC%', 'Transcription\nlevel']
    yter=1592477
    yori=3711828
    
    #Prepare GCSs data to bar-compatible format.
    GCSs_data_bared={}
    
    for ab, ab_ar in GSCs_data.items():
        bar_ar=[]
        for i in range(len(bins)-1):
            ab_num=0
            if bins[i]>=0:
                for gcs in ab_ar:
                    if bins[i+1]>gcs>=bins[i]:
                        ab_num+=1
                bar_ar.append(ab_num)
            else:
                for gcs in ab_ar: 
                    if bins[i+1]>gcs>=0 or genome_len>gcs>=bins[i]+genome_len:
                        ab_num+=1
                bar_ar.append(ab_num)
        bar_ar.append(bar_ar[0])
        print(ab, bar_ar)
        GCSs_data_bared[ab]=bar_ar
    
    #Compute confident intervals for GCSs number fall into MDs.
    MDs_confident_intervals={}  
    for ab, ab_ar in GSCs_data.items():
        upper_edge=[]
        lower_edge=[]
        print(MDs_lengths)
        for MD_len in MDs_lengths:
            upper_edge.append(binom.interval(0.999, len(ab_ar), MD_len/genome_len)[1])
            lower_edge.append(binom.interval(0.999, len(ab_ar), MD_len/genome_len)[0])
        upper_edge=[upper_edge[0]]+upper_edge+[upper_edge[0]]+[upper_edge[0]]
        lower_edge=[lower_edge[0]]+lower_edge+[lower_edge[0]]+[lower_edge[0]]
        MDs_confident_intervals[ab]=[lower_edge, upper_edge]
        print(ab, lower_edge)
        print(ab, upper_edge)
    
    #GCSs data plotting.
    fig, plot_names=plt.subplots(7,1,figsize=(11,15), dpi=100)
    print(bins)
    position=[0]+bins[1:]
    print(position)
    bin_width=[]
    position_bw=position+[genome_len]
    for j in range(len(position_bw)-1):
        bin_width.append(position_bw[j+1]-position_bw[j])
    print(bin_width)
    position_centre=[]
    for j in range(len(bin_width)):
        position_centre.append(int(position[j]+(bin_width[j]/2)))
    position_centre=[0]+position_centre+[genome_len]
    
    i=0
    Histo_comp_dict={} #Will contain computed histogramm data (bins and values)
    for key, value in GCSs_data_bared.items():
        plot_names[i].set_xlim(0, genome_len)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticks([yter, yori], minor=True)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        plot_names[i].locator_params(axis='y', nbins=6)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        plot_names[i].bar(position, value, bin_width, color=colors[i], linewidth=1, edgecolor='black', align='edge', zorder=1) #Barplot for GCSs number
        plot_names[i].plot(position_centre, MDs_confident_intervals[key][0], linestyle=":", color="black", linewidth=1, zorder=8) #Lower confidential border
        plot_names[i].plot(position_centre, MDs_confident_intervals[key][1], linestyle=":", color="black", linewidth=1, zorder=9) #Upper confidential border
        plot_names[i].plot(position_centre, MDs_confident_intervals[key][0], marker="_", color="black", linewidth=0, markersize=15, zorder=6) #Lower confidential border
        plot_names[i].plot(position_centre, MDs_confident_intervals[key][1],  marker="_", color="black", linewidth=0, markersize=15, zorder=7) #Upper confidential border
        plot_names[i].fill_between(position_centre, MDs_confident_intervals[key][0], MDs_confident_intervals[key][1], facecolor='grey', alpha=0.3, zorder=10) #Fill confident interval        
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(Y_labels[i], size=22, labelpad=8, rotation=90)
        i+=1
        
    #Score, GC, Transcription plotting.
    for key, value in Non_GCSs_data.items():
        value.append(value[0])
        plot_names[i].set_xlim(0, genome_len)
        if key=="GC":
            plot_names[i].set_ylim(45, max(value)+2)
        elif key=="Score":
            plot_names[i].set_ylim(min(value)-0.2, -1.5)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticks([yter, yori], minor=True)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        plot_names[i].locator_params(axis='y', nbins=6)
        plot_names[i].bar(position, value, bin_width, color=colors[i], linewidth=1, edgecolor='black', align='edge')
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(Y_labels[i], size=22, labelpad=8, rotation=90) 
        i+=1
    plt.tight_layout()
    fig.savefig(path_out+"GCSs_num_score_GC133_transcription_distrib_thr_genome.png", figsize=(11,15), dpi=400) 
    
    #GCSs data plotting for Cfx, Micro, and Oxo only.
    GSCs_data_main={'Cfx': GCSs_data_bared['Cfx'], 'Micro': GCSs_data_bared['Micro'], 'Oxo': GCSs_data_bared['Oxo']}
    Y_labels_main=['Cfx GCSs', 'Micro GCSs', 'Oxo GCSs']
    colors_main=['#7FCE79', '#ff878b', '#8991ff']
    fig, plot_names=plt.subplots(3,1,figsize=(11,7), dpi=100)
    i=0
    for key, value in GSCs_data_main.items():
        plot_names[i].set_xlim(0, genome_len)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticks([yter, yori], minor=True)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        plot_names[i].locator_params(axis='y', nbins=6)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        plot_names[i].bar(position, value, bin_width, color=colors_main[i], linewidth=1, edgecolor='black', align='edge')
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(Y_labels_main[i], size=22, labelpad=8, rotation=90)
        i+=1    
    plt.tight_layout()
    fig.savefig(path_out+"GCSs_number_Cfx_Micro_Oxo_distrib_thr_genome.png", figsize=(11,7), dpi=400) 
    return GCSs_data_bared, Non_GCSs_data


#######
#For evenly distributed bars.
#Computes following correlations: (Cfx GCSs data vs transcription), (Cfx GCSs data vs score), 
#(RifCfx GCSs data vs transcription), (RifCfx GCSs data vs score) 
#######

def track_corr(GCSs_histo_comp_dict, Non_GCSs_data):
    Cfx_GCSs=np.array(GCSs_histo_comp_dict['Cfx'][0]).tolist()
    RifCfx_GCSs=np.array(GCSs_histo_comp_dict['RifCfx'][0]).tolist()
    Micro_GCSs=np.array(GCSs_histo_comp_dict['Micro'][0]).tolist()
    Oxo_GCSs=np.array(GCSs_histo_comp_dict['Oxo'][0]).tolist()
    Transcription=Non_GCSs_data['Transcription']
    Score=Non_GCSs_data['Score']
    print('Paerson correlation (Cfx, transcription) for even bars: ' + str(pearsonr(Cfx_GCSs, Transcription)))
    print('Paerson correlation (RifCfx, transcription) for even bars: ' + str(pearsonr(RifCfx_GCSs, Transcription)))
    print('Paerson correlation (Cfx, score) for even bars: ' + str(pearsonr(Cfx_GCSs, Score)))
    print('Paerson correlation (RifCfx, score) for even bars: ' + str(pearsonr(RifCfx_GCSs, Score)))    
    return

#######
#For chromosomal macrodomains as bars.
#Computes following correlations: (Cfx GCSs data vs transcription), (Cfx GCSs data vs score), 
#(RifCfx GCSs data vs transcription), (RifCfx GCSs data vs score) 
#######

def track_corr_MDs(GCSs_data_bared, Non_GCSs_data_MDs):
    Transcription=Non_GCSs_data_MDs['Transcription']
    Score=Non_GCSs_data_MDs['Score']
    print('Paerson correlation (Cfx, transcription) for macrodomains: ' + str(pearsonr(GCSs_data_bared['Cfx'], Transcription)))
    print('Paerson correlation (RifCfx, transcription) for macrodomains: ' + str(pearsonr(GCSs_data_bared['RifCfx'], Transcription)))
    print('Paerson correlation (Cfx, score) for macrodomains: ' + str(pearsonr(GCSs_data_bared['Cfx'], Score)))
    print('Paerson correlation (RifCfx, score) for macrodomains: ' + str(pearsonr(GCSs_data_bared['RifCfx'], Score)))    
    return


#######
#Wrapps all the functions together.
#######

def wrap_the_functions(input_dict, broadpeak_path, score_path, GC_path, transcription_path, NAPs_input_dict, path_out, genome_len):
    #Even distribution of bins across the genome.
    GSCs_data=trusted_GCSs_parsing(input_dict)
    Non_GCSs_data_even=Prepare_non_GCSs_data(score_path, GC_path, transcription_path, np.linspace(0, genome_len, 11).tolist())
    GCSs_histo_comp_dict_even=Plot_the_distribution(GSCs_data, Non_GCSs_data_even, np.linspace(0, genome_len, 11).tolist(), genome_len, path_out+'Even_')
    track_corr(GCSs_histo_comp_dict_even, Non_GCSs_data_even)
    
    #Chromosomal macrodomains as genome bins.
    Macrodomains_descr=broadpeak_pars(broadpeak_path, genome_len)
    Non_GCSs_data_MDs=Prepare_non_GCSs_data(score_path, GC_path, transcription_path, Macrodomains_descr['Bins'])
    Bars=Plot_the_distribution_MDs(GSCs_data, Non_GCSs_data_MDs, Macrodomains_descr['Bins'], Macrodomains_descr['Lengths'], genome_len, path_out+'MDs_')
    track_corr_MDs(Bars[0], Bars[1])
    
    #Plot NAPs binding data.
    NAPs_data_even=BroadPeak_to_bins(NAPs_input_dict, np.linspace(0, genome_len, 11).tolist())[1]
    Plot_NAPs_distribution(NAPs_data_even, np.linspace(0, genome_len, 11).tolist(), path_out+"NAPs_binding_distrib_thr_genome_even.png", genome_len, 'even')
    NAPs_data_MDs=BroadPeak_to_bins(NAPs_input_dict, Macrodomains_descr['Bins'])[1]
    Plot_NAPs_distribution(NAPs_data_MDs, Macrodomains_descr['Bins'], path_out+"NAPs_binding_distrib_thr_genome_MDs.png", genome_len, 'MDs')
    return
    
wrap_the_functions(path_to_GCSs_files, Macrodomains_path, Score_data_path, GC_data_path, Transcription_data_path, NAPs_inpath, Plot_path_out, Genome_length)

print('Script ended its work succesfully!') 