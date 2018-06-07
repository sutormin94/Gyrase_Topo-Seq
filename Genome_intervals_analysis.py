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
path_to_TUs_sets={'All genes': '',
                  'High transcription level genes': '',
                  'Low transcription level genes': '',
                  'All operons': '',
                  'High transcription level operons': '',
                  'Low transcription level operons': '',
                  'Top long and active operons': '',
                  '16S operons': ''}
#Output for 16S operons analysis
Output_GCSs_num_16S=''

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
#GCSs association with rRNA operons.
#######

def rRNA_association(GCSs_sets_dict, operons_16S, window_width, path_out_16S):
    fileout=open(path_out_16S, 'w')
    fileout.write("Operon_ID\tGenes_ID\tStart\tEnd\tStrand\tExpression\tCondition\tGCSs in US\tGCSs in GB\tGCSs in DS\tSo on\n")    
    operons_stat_ds={} #{condition : [number of GCSs, mean N3E, mean score, len of regions]}
    ds_regions_len=0
    for j in operons_16S: #j - particular operon info (dictionary)  
        fileout.write(j['TUID'] + '\t' + j['TU name'] + '\t' + str(j['Start']) + '\t' + str(j['End']) + '\t' + j['Strand'] + '\t' + str(j['Transcription level']) + '\t')      
        ds_regions_len+=window_width
        for a, s in GCSs_sets_dict.items(): #a - Topo-Seq condition, s - corresponding set of GCSs.
            if a not in operons_stat_ds:
                operons_stat_ds[a]=[0, [], []] #[number of GCSs, mean N3E, mean score]
            fileout.write(a + '\t')
            stats=[0, 0, 0] #US GB DS
            N3E=[[], [], []] #US GB DS
            score=[[], [], []] #US GB DS
            for k, v in s.items(): #k - GCS coordinate, v - N3E and score.
                if j['Start']>k>j['Start']-window_width:
                    if j['Strand']=='-':
                        stats[2]+=1
                        N3E[2].append(v[0])
                        score[2].append(v[1])
                    elif j['Strand']=='+':
                        stats[0]+=1   
                        N3E[0].append(v[0])
                        score[0].append(v[1])
                elif j['End']+window_width>k>j['End']:
                    if j['Strand']=='-':
                        stats[0]+=1
                        N3E[0].append(v[0])
                        score[0].append(v[1])
                    elif j['Strand']=='+':
                        stats[2]+=1 
                        N3E[2].append(v[0])
                        score[2].append(v[1])
                elif j['End']>k>j['Start']:
                    stats[1]+=1
                    N3E[1].append(v[0])
                    score[1].append(v[1])
            fileout.write(str(stats[0]) + '\t' + str(stats[1]) + '\t' + str(stats[2]) + '\t')
            operons_stat_ds[a][0]+=stats[2]
            operons_stat_ds[a][1]+=N3E[2]
            operons_stat_ds[a][2]+=score[2]
        fileout.write('\n')
    for a, v in operons_stat_ds.items():
        v.append(ds_regions_len)
    fileout.close()
    return operons_stat_ds

#######
#Statistical analysis of intevals:
#1) Number of GCSs under/overrepresentation - binomail test.
#2) Interval GCSs mean N3E value vs all GCSs mean N3E - t-test.
#3) Interval GCSs mean score value vs all GCSs mean score - t-test.
#4) Interval mean score vs genome mean score - t-test.
#######

def interval_stat_analysis(GCSs_sets_dict, intervals_GCSs_dict, intervals, score_data, window_width, path_out):
    genome_len=4647454
    #Prepares lists of GCSs N3E and scores for statistical tests.
    GCSs_values_dict={}
    for a, s in GCSs_sets_dict.items():
        if a not in GCSs_values_dict:
            GCSs_values_dict[a]=[[s[0]], [s[1]]]
        else:
            GCSs_values_dict[a][0].append(s[0])
            GCSs_values_dict[a][1].append(s[1])
    fileout=open(path_out, 'r')
    #Performes t-test for comparison of intervals GCSs N3E and scores with all GCSs N3E and scores.
    #Performes binomail test for enrichment estimation of GCSs fall into intervals.
    fileout.write('Test\tAntibiotic\tValue\tInterval\tOverall\tp-value\n')
    for a, s in GCSs_sets_dict.items():
        N3E_stat=stats.ttest_ind(s[0], intervals_GCSs_dict[a][1]) #N3E
        score_stat=stats.ttest_ind(s[1], intervals_GCSs_dict[a][2]) #Score
        fileout.write('t-test\t' + a + '\tN3E\t' + str(np.mean(intervals_GCSs_dict[a][1])) + '\t' + 
                      str(np.mean(s[0])) + '\t' + str(N3E_stat[1]) + '\t' + 't-statistic: ' + str(N3E_stat[0]) + '\n')
        fileout.write('t-test\t' + a + '\tScore\t' + str(np.mean(intervals_GCSs_dict[a][2])) + '\t' + 
                      str(np.mean(s[1])) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(score_stat[0]) + '\n')
        GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][0], len(s[0]), intervals_GCSs_dict[a][3]/genome_len) #GSCs number
        fileout.write('binomial test\t' + a + '\tNumber of GCSs\t' + str(intervals_GCSs_dict[a][0]) + '\t' + 
                      str(len(s[0])) + '\t' + str(GCSs_number_stat) + '\n')
    #Performes t-test for comparison of intervals score with overall score.        
    intervals_values=[]
    for i in intervals:
        if i['Strand']=='+':
            intervals_values+=score_data[i['End']:i['End']+window_width]
        elif i['Strand']=='-':
            intervals_values+=score_data[i['Start']-window_width:i['Start']]
    intervals_score_stat=stats.ttest_ind(intervals_values, score_data) #Score
    fileout.write('t-test\t' + 'NA' + '\tInterval score\t' + str(np.mean(intervals_values)) + '\t' + 
                      str(np.mean(score_data)) + '\t' + str(intervals_score_stat[1]) + 't-statistic: ' + 
                      str(intervals_score_stat[0]) + '\n')
    fileout.close() 
    return

#######
#Add for:
#1) Other TUs sets
#2) BIMEs-2
#3) NAPs sites
#######



print('Script ended its work succesfully!') 