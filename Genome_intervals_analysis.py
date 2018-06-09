###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script analysis sets of genome intervals (transcription units, REPs, BIMEs-2, Fis sites, H-NS sites, MatP sites, etc.)
#for the enrichment of GCSs, compares their N3E and score with mean GCSs N3E and score, 
#compares intarvals mean score with genome mean score.
###############################################

#######
#Packages to be imported.
#######

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import scipy
from scipy.stats import binom
from scipy.stats.stats import pearsonr
from Bio import SeqIO
from Bio.SeqUtils import GC

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
Score_path=''
#Input data - sets of transcription units.
path_to_TUs_sets={'All_genes': '',
                  'High_transcription_level_genes': '',
                  'Low_transcription_level_genes': '',
                  'All_operons': '',
                  'High_transcription_level_operons': '',
                  'Low_transcription_level_operons': '',
                  'Top_long_and_active_operons': '',
                  '16S_operons': ''}
#Path for TU analysis output.
TU_analysis_outpath=''
#Input data - sets of intervals.
path_to_intervals_sets={'BIMEs1': '',
                  'BIMEs2': '',
                  'IHF sites': '',
                  'H-NS sites ': '',
                  'Fis sites': '',
                  'MatP sites': ''}
#Path for intervals analysis output.
Intervals_analysis_outpath=''

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


###############################################
#Group of functions for TUs analysis.
###############################################


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
#GCSs association with TUs.
#######

def TU_association(GCSs_sets_dict, TU_set, set_type, window_width, path_out, set_name):
    fileout=open(path_out + set_name + 'numbers_of_associated_GCSs.txt', 'w')
    if set_type=='16S operons':
        fileout.write("Operon_ID\tGenes_ID\tStart\tEnd\tStrand\tExpression\tCondition\tGCSs in US\tGCSs in GB\tGCSs in DS\tSo on\n") 
    elif set_type=='operons':
        fileout.write("Operon_ID\tGenes_ID\tStart\tEnd\tStrand\tExpression\tCondition\tGCSs in USUS\tGCSs in USGB\tGCSs in GBDS\tGCSs in DSDS\tSo on\n")  
    elif set_type=='genes':
        fileout.write("Gene_ID\tGene_name\tStart\tEnd\tStrand\tExpression\tCondition\tGCSs in USUS\tGCSs in USGB\tGCSs in GBDS\tGCSs in DSDS\tSo on\n")  
    else:
        print('Unknown type of the set! Possible types: 16S_operons, operons, genes.')
        return
    operons_stat_ds={} #{condition : [number of GCSs, mean N3E, mean score, len of regions]}
    ds_regions_len=0
    for j in TU_set: #j - particular operon info (dictionary)  
        fileout.write(j['TUID'] + '\t' + j['TU name'] + '\t' + str(j['Start']) + '\t' + str(j['End']) + '\t' + j['Strand'] + '\t' + str(j['Transcription level']) + '\t')      
        ds_regions_len+=window_width
        for a, s in GCSs_sets_dict.items(): #a - Topo-Seq condition, s - corresponding set of GCSs.
            if a not in operons_stat_ds:
                if set_type=='16S operons':
                    operons_stat_ds[a]=[0, [], []] #[number of GCSs, N3E values, score values]
                elif set_type=='operons' or set_type=='genes':
                    operons_stat_ds[a]=[0, [], [], 0, [], [], 0, [], [], 0, [], []] #[number of GCSs, N3E values, score values] for USUS, USGB, GBDS, DSDS correspondingly.
            fileout.write(a + '\t')
            stats=[0, 0, 0, 0, 0] #USUS USGB GB GBDS DSDS
            N3E=[[], [], [], [], []] #USUS USGB GB GBDS DSDS
            score=[[], [], [], [], []] #USUS USGB GB GBDS DSDS
            for k, v in s.items(): #k - GCS coordinate, v - N3E and score.
                if j['Start']>k>j['Start']-window_width and j['Strand']=='-': #DSDS
                    stats[4]+=1
                    N3E[4].append(v[0])
                    score[4].append(v[1])
                elif j['Start']>k>j['Start']-window_width and j['Strand']=='+': #USUS
                    stats[0]+=1   
                    N3E[0].append(v[0])
                    score[0].append(v[1])
                elif j['End']+window_width>k>j['End'] and j['Strand']=='-': #USUS
                    stats[0]+=1
                    N3E[0].append(v[0])
                    score[0].append(v[1])
                elif j['End']+window_width>k>j['End'] and j['Strand']=='+': #DSDS
                    stats[4]+=1 
                    N3E[4].append(v[0])
                    score[4].append(v[1])
                elif j['End']>k>j['Start']: #GB
                    stats[2]+=1
                    N3E[2].append(v[0])
                    score[2].append(v[1])
                elif j['Start']+window_width>k>j['Start'] and j['Strand']=='-': #GBDS
                    stats[3]+=1
                    N3E[3].append(v[0])
                    score[3].append(v[1])                    
                elif j['Start']+window_width>k>j['Start'] and j['Strand']=='+': #USGB
                    stats[1]+=1
                    N3E[1].append(v[0])
                    score[1].append(v[1])
                elif j['End']>k>j['End']-window_width and j['Strand']=='-': #USGB
                    stats[1]+=1
                    N3E[1].append(v[0])
                    score[1].append(v[1])  
                elif j['End']>k>j['End']-window_width and j['Strand']=='+': #GBDS
                    stats[3]+=1
                    N3E[3].append(v[0])
                    score[3].append(v[1]) 
            if set_type=='16S operons':
                fileout.write(str(stats[0]) + '\t' + str(stats[2]) + '\t' + str(stats[4]) + '\t') #US, GB, DS
                #DS
                operons_stat_ds[a][0]+=stats[4]
                operons_stat_ds[a][1]+=N3E[4]
                operons_stat_ds[a][2]+=score[4]                
            elif set_type=='operons' or set_type=='genes':
                fileout.write(str(stats[0]) + '\t' + str(stats[1]) + '\t' + str(stats[3]) + '\t' + str(stats[4]) + '\t') #USUS, USGB, GBDS, DSDS
                #USUS
                operons_stat_ds[a][0]+=stats[0]
                operons_stat_ds[a][1]+=N3E[0]
                operons_stat_ds[a][2]+=score[0] 
                #USGB
                operons_stat_ds[a][3]+=stats[1]
                operons_stat_ds[a][4]+=N3E[1]
                operons_stat_ds[a][5]+=score[1]                
                #GBDS
                operons_stat_ds[a][6]+=stats[3]
                operons_stat_ds[a][7]+=N3E[3]
                operons_stat_ds[a][8]+=score[3]                
                #DSDS
                operons_stat_ds[a][9]+=stats[4]
                operons_stat_ds[a][10]+=N3E[4]
                operons_stat_ds[a][11]+=score[4]                                
        fileout.write('\n')
    for a, v in operons_stat_ds.items():
        v.append(ds_regions_len)
    fileout.close()
    return operons_stat_ds #{condition : [##number of GCSs, N3E values, score values##*1-4, len of regions]}

#######
#Statistical analysis of intevals:
#1) Number of GCSs under/overrepresentation - binomail test.
#2) Interval GCSs mean N3E value vs all GCSs mean N3E - t-test.
#3) Interval GCSs mean score value vs all GCSs mean score - t-test.
#4) Interval mean score vs genome mean score - t-test.
#######

def TU_interval_stat_analysis(GCSs_sets_dict, intervals_GCSs_dict, intervals, score_data, window_width, set_type, path_out, set_name):
    genome_len=4647454
    #Prepares lists of GCSs N3E and scores for statistical tests.
    GCSs_values_dict={}
    for a, s in GCSs_sets_dict.items():
        for k, v in s.items():
            if a not in GCSs_values_dict:
                GCSs_values_dict[a]=[[v[0]], [v[1]]]
            else:
                GCSs_values_dict[a][0].append(v[0])
                GCSs_values_dict[a][1].append(v[1])
            
    fileout=open(path_out + set_name + 'compartments_statistics_GCSs_number_N3E_score.txt', 'w')
    #Performes t-test for comparison of intervals GCSs N3E and scores with all GCSs N3E and scores.
    #Performes binomail test for enrichment estimation of GCSs fall into intervals.
    fileout.write('Test\tAntibiotic\t')
    for i in range((len(next(iter(intervals_GCSs_dict.values())))-1)/3):
        fileout.write('Examined value\tValue (intervals)\tValue (overall)\tp-value\tAdditional\t')
    fileout.write('\n') 
    
    if set_type=='16S operons':
        for a, ns in GCSs_values_dict.items():
            N3E_stat=stats.ttest_ind(ns[0], intervals_GCSs_dict[a][1]) #N3E DS
            fileout.write('t-test\t' + a + '\tDS N3E\t' + str(np.mean(intervals_GCSs_dict[a][1])) + '\t' + 
                          str(np.mean(ns[0])) + '\t' + str(N3E_stat[1]) + '\t' + 't-statistic: ' + str(N3E_stat[0]) + '\n')            
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][2]) #Score DS
            fileout.write('t-test\t' + a + '\tDS Score\t' + str(np.mean(intervals_GCSs_dict[a][2])) + '\t' + 
                          str(np.mean(ns[1])) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(score_stat[0]) + '\n')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][0], len(ns[0]), intervals_GCSs_dict[a][3]/genome_len) #GSCs number DS
            fileout.write('binomial test\t' + a + '\tDS Number of GCSs\t' + str(intervals_GCSs_dict[a][0]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\n')
    
    elif set_type=='operons' or set_type=='genes':
        for a, ns in GCSs_values_dict.items():
            #N3E
            N3E_stat=stats.ttest_ind(ns[0], intervals_GCSs_dict[a][1]) #N3E USUS
            fileout.write('t-test\t' + a + '\tUSUS N3E\t' + str(np.mean(intervals_GCSs_dict[a][1])) + '\t' + 
                          str(np.mean(ns[0])) + '\t' + str(N3E_stat[1]) + '\t' + 't-statistic: ' + str(N3E_stat[0]) + '\t')   
            N3E_stat=stats.ttest_ind(ns[0], intervals_GCSs_dict[a][4]) #N3E USGB
            fileout.write('USGB N3E\t' + str(np.mean(intervals_GCSs_dict[a][4])) + '\t' + 
                          str(np.mean(ns[0])) + '\t' + str(N3E_stat[1]) + '\t' + 't-statistic: ' + str(N3E_stat[0]) + '\t')               
            N3E_stat=stats.ttest_ind(ns[0], intervals_GCSs_dict[a][7]) #N3E GBDS
            fileout.write('GBDS N3E\t' + str(np.mean(intervals_GCSs_dict[a][7])) + '\t' + 
                          str(np.mean(ns[0])) + '\t' + str(N3E_stat[1]) + '\t' + 't-statistic: ' + str(N3E_stat[0]) + '\t')   
            N3E_stat=stats.ttest_ind(ns[0], intervals_GCSs_dict[a][10]) #N3E DSDS
            fileout.write('DSDS N3E\t' + str(np.mean(intervals_GCSs_dict[a][10])) + '\t' + 
                          str(np.mean(ns[0])) + '\t' + str(N3E_stat[1]) + '\t' + 't-statistic: ' + str(N3E_stat[0]) + '\n')               
            #Score
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][2]) #Score USUS
            fileout.write('t-test\t' + a + '\tUSUS Score\t' + str(np.mean(intervals_GCSs_dict[a][2])) + '\t' + 
                          str(np.mean(ns[1])) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(score_stat[0]) + '\t')
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][5]) #Score USGB
            fileout.write('USGB Score\t' + str(np.mean(intervals_GCSs_dict[a][5])) + '\t' + 
                          str(np.mean(ns[1])) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(score_stat[0]) + '\t')
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][8]) #Score GBDS
            fileout.write('GBDS Score\t' + str(np.mean(intervals_GCSs_dict[a][8])) + '\t' + 
                          str(np.mean(ns[1])) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(score_stat[0]) + '\t')
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][11]) #Score DSDS
            fileout.write('DSDS Score\t' + str(np.mean(intervals_GCSs_dict[a][11])) + '\t' + 
                          str(np.mean(ns[1])) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(score_stat[0]) + '\n')            
            #GCSs number
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][0], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len) #GSCs number USUS
            fileout.write('binomial test\t' + a + '\tUSUS Number of GCSs\t' + str(intervals_GCSs_dict[a][0]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\tNA\t')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][3], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len) #GSCs number USGB
            fileout.write('USGB Number of GCSs\t' + str(intervals_GCSs_dict[a][3]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\tNA\t')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][6], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len) #GSCs number GBDS
            fileout.write('GBDS Number of GCSs\t' + str(intervals_GCSs_dict[a][6]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\tNA\t')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][9], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len) #GSCs number DSDS
            fileout.write('DSDS Number of GCSs\t' + str(intervals_GCSs_dict[a][9]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\tNA\n')            
    fileout.write('\n')
    
    #Performes t-test for comparison of intervals score with overall score.
    if set_type=='16S operons': #DS data only
        intervals_values=[]
        for i in intervals:
            if i['Strand']=='+':
                intervals_values+=score_data[i['End']:i['End']+window_width]
            elif i['Strand']=='-':
                intervals_values+=score_data[i['Start']-window_width:i['Start']]
        intervals_score_stat=stats.ttest_ind(intervals_values, score_data) #Score
        fileout.write('t-test\t' + 'NA' + '\t DS Interval score\t' + str(np.mean(intervals_values)) + '\t' + 
                      str(np.mean(score_data)) + '\t' + str(intervals_score_stat[1]) + 't-statistic: ' + 
                      str(intervals_score_stat[0]) + '\n')
    
    elif set_type=='operons' or set_type=='genes':
        intervals_values=[[], [], [], []] #USUS, USGB, GBDS, DSDS
        for i in intervals:
            if i['Strand']=='+':
                intervals_values[0]+=score_data[i['Start']-window_width:i['Start']] #USUS
                intervals_values[1]+=score_data[i['Start']:i['Start']+window_width] #USGB
                intervals_values[2]+=score_data[i['End']-window_width:i['End']] #GBDS
                intervals_values[3]+=score_data[i['End']:i['End']+window_width] #DSDS
            elif i['Strand']=='-':
                intervals_values[0]+=score_data[i['Start']-window_width:i['Start']] #DSDS
                intervals_values[1]+=score_data[i['Start']:i['Start']+window_width] #GBDS
                intervals_values[2]+=score_data[i['End']-window_width:i['End']] #USGB
                intervals_values[3]+=score_data[i['End']:i['End']+window_width] #USUS
        intervals_score_stat=[]
        for i in intervals_values:
            intervals_score_stat.append(stats.ttest_ind(i, score_data))
        fileout.write('t-test\t' + 'NA' + 
                      '\t USUS Interval score\t' + str(np.mean(intervals_values[0])) + '\t' + str(np.mean(score_data)) + '\t' + 
                      str(intervals_score_stat[0][1]) + 't-statistic: ' + str(intervals_score_stat[0][0]) + '\t' + 
                      '\t USGB Interval score\t' + str(np.mean(intervals_values[1])) + '\t' + str(np.mean(score_data)) + '\t' + 
                      str(intervals_score_stat[1][1]) + 't-statistic: ' + str(intervals_score_stat[1][0]) + '\t' + 
                      '\t GBDS Interval score\t' + str(np.mean(intervals_values[2])) + '\t' + str(np.mean(score_data)) + '\t' + 
                      str(intervals_score_stat[2][1]) + 't-statistic: ' + str(intervals_score_stat[2][0]) + '\t' + 
                      '\t DSDS Interval score\t' + str(np.mean(intervals_values[3])) + '\t' + str(np.mean(score_data)) + '\t' + 
                      str(intervals_score_stat[3][1]) + 't-statistic: ' + str(intervals_score_stat[3][0]) + '\n')      
    fileout.close() 
    return

#######
#Number of GCSs statistical analysis for those associated with TUs, additional normalization step is added. 
#parameters_example=[Number of GCSs for particular condition, Total number of TUs, Number of TUs in the partial set]
#######

def GCSs_num_normalization(all_TUs, part_TU, parameters):
    GCSs_set_exp_interval_dict={} 
    for a, s in all_TUs.items():
        GCSs_all=all_TUs[a][0]+all_TUs[a][3]+all_TUs[a][6]+all_TUs[a][9] #Total number of GCSs observed to be associated with the full set of TUs.
        GCSs_all_exp_compartment=float(GCSs_all)/4 #Expected number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for full set of TUs.
        GCSs_set_exp_interval=GCSs_all_exp_compartment*parameters[2]/parameters[1] #Expected number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for partial set of TUs.
        GCSs_set_exp_interval_dict[a]=[GCSs_set_exp_interval] 
        for j in range(((len(s))-1)/3):
            GCSs_set_exp_interval_dict[a].append(part_TU[a][j*3]) #Number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for partial set of TUs.
            GCSs_set_exp_interval_dict[a].append(binom.cdf(part_TU[a][j*3], GCSs_set_exp_interval*4, 0.25)) #p-value of binomial test for the number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for partial set of TUs.
            GCSs_set_exp_interval_dict[a].append(float(part_TU[a][j*3])*100000/(parameters[0]*parameters[2])) #Normalized number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for partial set of TUs.
    return GCSs_set_exp_interval_dict #[GCSs expected] + [GCSs obs, p-value, GCSs norm]*[USUS, USGB, GBDS, DSDS]

#######
#Writes GCSs numbers information to file: Number of GCSs expected, Number of GCSs observed, p-value - binomial test, Number of GCSs normalized.
#######

def write_GCSs_norm(GCSs_set_exp_interval_dict, path_out, set_name):
    fileout=open(path_out + set_name + 'normalized_GCSs_numbers_and_statistics.txt', 'w')
    fileout.write('Condition\tCompartment\tNumber of GCSs expected\tNumber of GCSs observed\tp-value\tNumber of GCSs normalized\n')
    Compartment_names=['USUS', 'USGB', 'GBDS', 'DSDS']
    for a, s in GCSs_set_exp_interval_dict.items():
        for i in range(len(Compartment_names)):
            fileout.write(a + '\t' + Compartment_names[i] + '\t' + str(s[0]) + '\t' + str(s[(i*3)+1]) + '\t' + str(s[(i*3)+2]) + '\t' + str(s[(i*3)+3]) +'\n')
    fileout.close()
    return

###############################################
#Group of functions for intervals analysis.
###############################################

#######
#BroadPeak-formalized intervals parsing (NAPs sites or BIMEs, etc) and filtering intervals that are not deleted.
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

def GCSs_in_intervals(GCSs_sets_dict,  intervals, path_out):
    fileout=open(path_out+'GCSs_associated_with_intervals_statistics.txt', 'w')
    fileout.write('Interval type\tCondition\tNumber of GCSs expected\tNumber of GCSs observed\tp-value\t' +
                  'GCSs N3E (intervals)\tGCSs N3E (overall)\tp-value\tt-statistic\t' + 
                  'GCSs score (intervals)\tGCSs score (overall)\tp-value\tt-statistic\n')
    #Correcting genome length.
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    del_len=0
    for i in deletions:
        del_len+=i[1]-i[0]
    genome_len=4647454
    genome_len_dc=genome_len-del_len  
    
    #Prepares lists of GCSs N3E and scores for statistical tests.
    GCSs_values_dict={}
    for a, s in GCSs_sets_dict.items():
        for k, v in s.items():
            if a not in GCSs_values_dict:
                GCSs_values_dict[a]=[[v[0]], [v[1]]]
            else:
                GCSs_values_dict[a][0].append(v[0])
                GCSs_values_dict[a][1].append(v[1])    
    
    #Finds GCSs associated with intervals.
    GCSs_associated_info={} #Dictionary contains element corresponds to interval sets.
    for k, v in intervals.items(): #Iterate different interval types.
        GCSs_associated_info[k]={}
        intervals_len=0
        for i in v: #Iterate intevals.
            intervals_len+=i[1]-i[0]      
            for a, s in GCSs_sets_dict.items(): #Iterate GCSs sets.
                if a not in GCSs_associated_info[k]:
                    GCSs_associated_info[k][a]=[0, [], []] #Dictionary contains elements corresponds to GCSs sets. [Number of GCSs, N3E values, Score values]
                for gcs, info in s.items(): #Iterate GCSs.
                    if i[1]>gcs>i[0]:
                        GCSs_associated_info[k][a][0]+=1
                        GCSs_associated_info[k][a][1].append(info[0])
                        GCSs_associated_info[k][a][2].append(info[1])
        
        for a, s in GCSs_associated_info[k].items():
            relative_intervals_len=float(intervals_len)/genome_len
            total_number_of_GCSs=len(GCSs_sets_dict[a])
            expected_number_of_GCSs_in_interval=relative_intervals_len*total_number_of_GCSs 
            GCSs_num_stat=binom.cdf(s[0], total_number_of_GCSs, relative_intervals_len) #Number of GCSs stat.
            GCSs_N3E_stat=stats.ttest_ind(s[1], GCSs_values_dict[a][0]) #GCSs N3E stat.
            GCSs_score_stat=stats.ttest_ind(s[2], GCSs_values_dict[a][1]) #GCSs score stat.
            fileout.write(k + '\t' + a + '\t' + str(expected_number_of_GCSs_in_interval) + '\t' + str(s[0]) + '\t' + str(GCSs_num_stat) + '\t' + 
                          str(np.mean(s[1])) + '\t' + str(np.mean(GCSs_values_dict[a][0])) + '\t' + str(GCSs_N3E_stat[1]) + '\t' + str(GCSs_N3E_stat[0]) + '\t' + 
                          str(np.mean(s[2])) + '\t' + str(np.mean(GCSs_values_dict[a][1])) + '\t' + str(GCSs_score_stat[1]) + '\t' + str(GCSs_score_stat[0]) + '\n')
    
    fileout.close()
    return GCSs_associated_info


print('Script ended its work succesfully!') 