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

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.stats import binom

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Input data - GCSs, TAB.
path_to_GCSs_files={'Cfx': '',
                    'RifCfx': '',
                    'Micro': '',
                    'Oxo': ''}

#Input data - score, WIG.
Score_data_path=''

#Input data - GC, WIG.
GC_data_path=''

#Input data - expression, TAB.
Transcription_data_path=''

#Output data: plot.
Plot_path_out=''

#######
#GCSs data parsing.
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
        for j in range(int(line[3])-int(line[2])):
            transcription[int(line[2])+j]=float(line[5])
    print('Whole genome average transcription: ' + str(sum(transcription)/len(transcription)))
    return transcription

#######
#Convert score/transcription/GC data to histogram.
#######

def bar_convert(ar):
    bar_ar=[]
    bin_width=int(4647454/9)
    for i in range(9):
        bar_ar.append(sum(ar[i*bin_width:(i+1)*bin_width])/bin_width)
    return bar_ar

#######
#Wrapps all the nonGCSs data together and prepare it for plotting.
#######

def Prepare_non_GCSs_data(score_path, GC_path, transcription_path):
    Data_dict={}
    Score_data=score_data_parser(score_path, 'Score')
    GC_data=score_data_parser(GC_path, 'GC')
    Trancription_data=transcription_data_parser(transcription_path)
    Data_dict['Score']=bar_convert(Score_data)
    Data_dict['GC']=bar_convert(GC_data)
    Data_dict['Transcription']=bar_convert(Trancription_data)
    return Data_dict

#######
#Plotting: distribution throughout the genome.
#######

def Plot_the_distribution(GSCs_data, Non_GCSs_data, path_out):
    #Parameters
    ticks1=[0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000]
    xticknames1=['', '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500']
    colors=['#ff545a', '#fa585a', '#49ba41', '#5762ff', '#ac5eff', '#50b3ff', '#ffd75e']
    plot_names=['plot1', 'plot2', 'plot3', 'plot4', 'plot5', 'plot6', 'plot7']
    Y_labels=['Cfx GCSs', 'RifCfx GCSs', 'Micro GCSs', 'Oxo GCSs', 'Score', 'GC%', 'Transcription\nlevel']
    yter=1592477
    yori=3711828
    #GCSs data plotting.
    fig, plot_names=plt.subplots(7,1,figsize=(11,15))
    i=0
    for key, value in GSCs_data.items():
        bins=np.linspace(0, 4647454, 9)
        plot_names[i].set_xlim(0, 4647454)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        conf_interval=stat(binom.interval(0.95, len(value), 1/9)[0], binom.interval(0.95, len(value), 1/9)[1])
        plot_names[i].set_yticks(conf_interval, minor=True)
        plot_names[i].yaxis.grid(True, which='minor', linewidth=0.4, linestyle='--', color='black')
        plot_names[i].fill_between(bins, conf_interval[0], conf_interval[1], facecolor='grey', alpha=0.3)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        plot_names[i].hist(value, bins, facecolor=colors[i], alpha=0.7, linewidth=1, edgecolor='black')
        plot_names[i].annotate('Origin \nof replication', xytext=(0, 50), textcoords='offset points', xy=(yori, 0), xycoords='data', color='black', arrowprops=dict(facecolor='black', shrink=0, width=0.3), size=16.5, weight="bold")
        plot_names[i].annotate('Terminators \n area', xytext=(0, 50), textcoords='offset points', xy=(yter, 0), xycoords='data', color='black', arrowprops=dict(facecolor='black', shrink=0, width=0.3), size=16.5, weight="bold")
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(Y_labels[i], size=22, labelpad=8, rotation=90)
        i+=1
    #Score, GC, Transcription plotting.
    bin_width=int(4647454/9)
    position=[]
    for i in range(9):
        position.append((i)*bin_width)   
    for key, value in Non_GCSs_data.items():
        bins=np.linspace(0, 4647454, 9)
        plot_names[i].set_xlim(0, 4647454)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        conf_interval=stat(binom.interval(0.95, len(value), 1/9)[0], binom.interval(0.95, len(value), 1/9)[1])
        plot_names[i].set_yticks(conf_interval, minor=True)
        plot_names[i].yaxis.grid(True, which='minor', linewidth=0.4, linestyle='--', color='black')
        plot_names[i].fill_between(bins, conf_interval[0], conf_interval[1], facecolor='grey', alpha=0.3)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        plot_names[i].bar(position, value, bin_width, color=colors[i], linewidth=1, edgecolor='black', align='edge')
        plot_names[i].annotate('Origin \nof replication', xytext=(0, 50), textcoords='offset points', xy=(yori, 0), xycoords='data', color='black', arrowprops=dict(facecolor='black', shrink=0, width=0.3), size=16.5, weight="bold")
        plot_names[i].annotate('Terminators \n area', xytext=(0, 50), textcoords='offset points', xy=(yter, 0), xycoords='data', color='black', arrowprops=dict(facecolor='black', shrink=0, width=0.3), size=16.5, weight="bold")
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(Y_labels[i], size=22, labelpad=8, rotation=90) 
        i+=1
    plt.tight_layout(pad=1.08, h_pad=2, w_pad=None, rect=None)
    plt.show()
    fig.savefig(path_out, dpi=320) 
    return

#######
#Wrapps all the functions together.
#######

def wrap_the_functions(input_dict, score_path, GC_path, transcription_path, path_out):
    GSCs_data=trusted_GCSs_parsing(input_dict)
    Non_GCSs_data=Prepare_non_GCSs_data(score_path, GC_path, transcription_path)
    Plot_the_distribution(GSCs_data, Non_GCSs_data, path_out)
    return
    
wrap_the_functions(path_to_GCSs_files, Score_data_path, GC_data_path, Transcription_data_path, Plot_path_out)

print('Script ended its work succesfully!') 