###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Script compares data from Cfx and RifCfx (conditions with transcription inhibited with rifampicin) experiments. It identifies GCSs
#shared between datasets, computes whether the signal (N3E) goes up or down as a response for transcription inhibition. Also it 
#analysis shared GCSs that fall into BIMEs or DS regions of rRNA operons to be associated with signal increase or decrease.
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
                    'RifCfx': ''}
#Input data - sets of transcription units.
path_to_TUs_sets={'16S_operons': ''}
#Input data - sets of intervals.
path_to_intervals_sets={'BIMEs1': '',
                  'BIMEs2': ''}
#Outpath for plot.
Plot_path_out=''


#######
#Trusted GCSs data parsing.
#######

def trusted_GCSs_parsing(input_dict):
    GCSs_sets_dict={}
    for k, v in input_dict.items():
        GCSs_dict={}
        filein=open(v, 'r')
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in ['GCSs_coordinate']:
                GCSs_dict[int(line[0])]=[float(line[1]), float(line[2])]
            else:
                continue
        filein.close()
        GCSs_sets_dict[k]=GCSs_dict
        print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(GCSs_dict)))
    return GCSs_sets_dict


#######
#Prepare lists of GCSs heights for Cfx and Cfx_Rif data. Write down values that higher for RifCfx.
#######


def find_shared_gcss(GCSs_sets_dict):
    #Shared GCSs identification
    shared_GCSs_array=[]
    for k, v in GCSs_sets_dict['Cfx'].items():
            if k in GCSs_sets_dict['RifCfx']:
                shared_GCSs_array.append([k, GCSs_sets_dict['RifCfx'][k][0]/v[0], GCSs_sets_dict['RifCfx'][k][1]]) #[GCSs coordinate, N3E ratio (RifCfx/Cfx), GCSs score]
    #GCSs sorting by N3E ratio
    shared_GCSs_array.sort(key=lambda tup: tup[1])
    #Prepare data for plotting
    j=0
    for i in shared_GCSs_array:
        i.append(j) #Serial number
        j+=1
    return shared_GCSs_array #array consist of elements: [GCSs coordinate, N3E ratio (RifCfx/Cfx), GCSs score, Serial number]


#######
#Parsing TUs sets.
#######

def TUs_parser(TUs_sets_path):
    TUs_sets={}
    for k, v in TUs_sets_path.items():
        filein=open(v, 'r')
        plus=0
        minus=0
        ar=[]
        for line in filein:
            tu={}
            line=line.rstrip().split('\t')
            if line[0] not in ["GeneID", "OperonID"]:
                tu['TUID']=str(line[0]) #TUID
                tu['TU name']=str(line[1]) #TU name/composition
                tu['Start']=int(line[2]) #Start
                tu['End']=int(line[3]) #End
                tu['Strand']=str(line[4]) #Strand +/-
                tu['Transcription level']=float(line[5]) #Transcription level
                tu['Description']=str(line[6]) #Description
                ar.append(tu)
                if line[4]=='+':
                    plus+=1
                elif line[4]=="-":
                    minus+=1
        filein.close()
        TUs_sets[k]=ar
        print('Number of TUs in forward for ' + str(k) + 'set: ' + str(plus))
        print('Number of TUs in reverse for ' + str(k) + 'set: ' + str(minus))
    return TUs_sets


#######
#Find GCSs associated with TUs DS region.
#######

def TU_association(GCSs_set, TU_set, window_width):
    #Find GCSs fall into DS compartments of TUs
    GCSs_in_DS=[] #GCSs fall into DS compartment
    for j in TU_set: #Iterate TUs 
        for k in GCSs_set: #Iterate GCSs
            if j['Start']>k[0]>j['Start']-window_width and j['Strand']=='-': #DSDS
                GCSs_in_DS.append(k)
            elif j['End']+window_width>k[0]>j['End'] and j['Strand']=='+': #DSDS
                GCSs_in_DS.append(k)
    #GCSs sorting by serial number ratio
    GCSs_in_DS.sort(key=lambda tup: tup[3])
    return GCSs_in_DS #array of elements: [GCSs coordinate, N3E ratio (RifCfx/Cfx), GCSs score, Serial number]


#######
#BroadPeak-formalized intervals parsing (BIMEs or NAPs sites, etc.) and filtering intervals that are not deleted.
#######

def broadpeak_pars(intervals_sets_path):
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    intervals_sets_dict={}
    for k, v in intervals_sets_path.items():
        filein=open(v, 'r')
        ar=[]
        for line in filein:
            line=line.rstrip().split('\t')
            int_start=int(line[1])
            int_end=int(line[2])
            del_check=0
            for j in range(len(deletions)):
                if deletions[j][1]>int_start>deletions[j][0] or deletions[j][1]>int_end>deletions[j][0]:
                    del_check=1
            if del_check==0:
                ar.append([int_start, int_end])            
        intervals_sets_dict[k]=ar
        print("Number of " + str(k) + "regions: " + str(len(ar)))
        filein.close()
    return intervals_sets_dict

#######
#GCSs in intervals.
#######   

def intervals_association(GCSs_set, intervals):   
    #Correcting genome length.
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    del_len=0
    for i in deletions:
        del_len+=i[1]-i[0]
    genome_len=4647454
    genome_len_dc=genome_len-del_len     
    
    #Finds GCSs associated with intervals
    GCSs_in_intervals=[]
    for i in intervals: #Iterate intervals
        for k in GCSs_set: #Iterate GCSs
            if i[1]>k[0]>i[0]:
                GCSs_in_intervals.append(k)
    #GCSs sorting by serial number ratio
    GCSs_in_intervals.sort(key=lambda tup: tup[3])
    return GCSs_in_intervals #array of elements: [GCSs coordinate, N3E ratio (RifCfx/Cfx), GCSs score, Serial number]

#######
#Calculates number of GCSs that increase and decrease signal when transcription is inhibited.
#######   

def comp_with_one(ar):
    increase=0
    decrease=0
    for i in ar:
        if i[1]>1:
            increase+=1
        else:
            decrease+=1
    return [increase, decrease]

####### 
#Performs binomial test to find preferential behaviour (to increase or to decrease).
####### 

def binom_test_preferent(all_GCSs, set_GCSs, set_type):
    indec_all=comp_with_one(all_GCSs)
    indec_set=comp_with_one(set_GCSs)
    print("Number of GCSs increase N3E when transcription is inhibited: " + str(indec_all[0]) + " out of " + str(indec_all[0]+indec_all[1]))
    prob_to_increase=indec_all[0]/(float(indec_all[0]+indec_all[1]))
    print("Probability a GCSs will increase N3E: " + str(prob_to_increase) + "\n")
    print("Number of GCSs fall into " + str(set_type) + " that increase N3E when transcription is inhibited: " + str(indec_set[0]) + " out of " + str(indec_set[0]+indec_set[1]))
    prob_of_distrib=binom.cdf(indec_set[0], indec_set[0]+indec_set[1], prob_to_increase)
    if prob_of_distrib>0.5:
        prob_of_distrib=1-prob_of_distrib
    print("Probability of this distribution: " + str(prob_of_distrib) + "\n")

####### 
#Delete GCSs repeated shared between several GCSs sets (left this GCSs in particular only set).
####### 

def exclude_repeats(ar1, ar2):
    elements_to_keep=[]
    for i in ar1: #Iterate set one
        check=0
        for j in ar2: #Itarate set two (which to exclude)
            if i[0]==j[0]:
                check=1
        if check==0:
            elements_to_keep.append(i)
    return elements_to_keep


#######
#Sorted N3E ratio (RifCfx/Cfx) curve plotting.
#######

def plot_N3E_ratio(dict_of_sets, pathout):
    fig = plt.figure(figsize=(16,6))                                                               
    plt1 = fig.add_subplot(1,1,1)  
    ticks1=[1]
    ticks2=[int(len(dict_of_sets['all_GCSs'][0])/2)]
    plt1.set_xlim(0, len(dict_of_sets['all_GCSs'][0]))
    N3E_ratio=[]
    for i in dict_of_sets['all_GCSs'][0]:
        N3E_ratio.append(i[1])
    plt1.set_ylim(0, max(N3E_ratio)+1)
    plt1.set_yticks(ticks1, minor=True)
    plt1.set_xticks(ticks2, minor=True)
    plt1.grid(True, axis='y', which='minor', linewidth=2, linestyle='--', alpha=0.9)
    plt1.grid(True, axis='x', which='minor', linewidth=2, linestyle='--', alpha=0.9)
    for k, v in dict_of_sets.items(): #Iterate GCSs sets.
        if k not in ['all_GCSs']: 
            for gcs in v[0]: #Iterate GCSs
                plt1.scatter(v[0][3], v[0][1], c=v[1][0], alpha=v[1][1], linewidths=v[1][2], s=v[1][3])
    
    plt1.scatter(620, 9.5, c='#74ffa0', alpha=1, linewidths=1, s=200)
    plt1.annotate('GCSs in BIMEs2', xytext=(635, 9.3), xy=(0, 0), xycoords='data', color='black', weight="bold", size=27)  
    plt1.scatter(620, 8.8, c='#ff8074', alpha=1, linewidths=1, s=200)
    plt1.annotate('GCSs in rRNA operons DS', xytext=(635, 8.6), xy=(0, 0), xycoords='data', color='black', weight="bold", size=27)  
    plt1.scatter(620, 8.1, c='#8b74ff', alpha=0.2, linewidths=0, s=100)
    plt1.annotate('Other sites', xytext=(635, 7.9), xy=(0, 0), xycoords='data', color='black', weight="bold", size=27)  
    plt1.annotate('Median', xytext=(535, 4.5), xy=(0, 0), xycoords='data', rotation=90, color='black', size=25)  

    xticknames1=[0, 100, 200, 300, 400, 500,600,700,800,900,1000,1100]
    yticknames1=[0,1,2,3,4,5,6,7,8,9,10,11]
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_xticks(xticknames1, minor=False)
    plt1.set_xlabel('Sites', size=33, labelpad=8, rotation=0)
    plt1.set_xticklabels(xticknames1)
    plt.setp(plt1.set_xticklabels(xticknames1), rotation=0, fontsize=25)
    plt1.set_yticklabels(yticknames1)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=25)
    plt1.set_ylabel('RifCfx/Cfx', size=33, labelpad=8, rotation=90)
    fig.savefig(pathout, dpi=320) 
    return

#######
#Wrapps all the functions together.
#######

def functions_wrapper(GCSs_dict_input, TUs_sets_path, intervals_sets_path, plot_output):
    #Parsing GCSs
    GCSs_sets_dict=trusted_GCSs_parsing(GCSs_dict_input)
    #Find GCSs shared between Cfx and RifCfx
    shared_GCSs_set=find_shared_gcss(GCSs_sets_dict)
    #Parsing TUs
    TU_set=TUs_parser(TUs_sets_path)['16S_operons']
    #Find shared GCSs associated with DS compartment of TUs
    window_width=5000
    GCSs_in_DS=TU_association(shared_GCSs_set, TU_set, window_width)
    #Parsing intervals
    intervals_set=broadpeak_pars(intervals_sets_path)['BIMEs2']
    #Find shared GCSs associated with intervals
    GCSs_in_intervals=intervals_association(shared_GCSs_set, intervals_set)
    #Revealing increase/decrease trends within a GCSs set
    binom_test_preferent(shared_GCSs_set, GCSs_in_DS, '16S operons DS')
    binom_test_preferent(shared_GCSs_set, GCSs_in_intervals, 'BIMEs2')
    #Filter GCSs sets
    GCSs_in_DS_only=exclude_repeats(GCSs_in_DS, GCSs_in_intervals)
    shared_GCSs_minus_intervals_DS=exclude_repeats(shared_GCSs_set, GCSs_in_intervals+GCSs_in_DS_only)
    #Prepare data for plotting
    dict_of_sets={'all_GCSs': [shared_GCSs_set, ['#8b74ff', 0.15, 0, 100]], 
                  'rRNA_DS': [GCSs_in_DS_only, ['#ff8074', 1, 1, 200]],
                  'BIMEs2': [GCSs_in_intervals, ['#74ffa0', 1, 1, 200]], 
                  'other': [shared_GCSs_minus_intervals_DS, ['#8b74ff', 0.15, 0, 100]]}
    plot_N3E_ratio(dict_of_sets, plot_output)
    return
    
print('Script ended its work succesfully!') 