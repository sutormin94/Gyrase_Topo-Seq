###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script analysis sets of genome intervals (transcription units - TUs, BIMEs-1, BIMEs-2, IHF sites, Fis sites, H-NS sites, MatP sites, etc.)
#for the enrichment of GCSs (binomial test), compares their N3E and score with mean GCSs N3E and score (t-test), 
#compares intervals mean score with genome mean score (t-test).
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
from scipy.stats import binom

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
Score_path="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\E_coli_w3110_G_Mu_score.wig"
#Input data - sets of transcription units.
path_to_TUs_sets={'All_genes': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_genes_expression.txt",
                  'High_transcription_level_genes': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_high_expression_genes_370.txt",
                  'Low_transcription_level_genes': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_low_expression_genes_370.txt",
                  'All_operons': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_operons_expression.txt",
                  'High_transcription_level_operons': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_high_expression_operons_186.txt",
                  'Low_transcription_level_operons': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_low_expression_operons_186.txt",
                  'Top_long_and_active_operons': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Top_long_and_active_DOOR_Mu_operons_expression.txt",
                  '16S_operons': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\DOOR_Mu_del_cor_16S_rRNA_operons.txt"}
#Length of the genome, bp
Genome_length=4647454
#Path for TU analysis output.
TU_analysis_outpath="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_association_with_TUs\\"
if not os.path.exists(TU_analysis_outpath):
    os.makedirs(TU_analysis_outpath)
#Input data - sets of intervals.
'''
path_to_intervals_sets={'BIMEs1': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\BIMEs1_coordinates.broadPeak",
                  'BIMEs2': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\BIMEs2_coordinates.broadPeak",
                  'IHFA': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\IHFA_sites_mid_exp_Prieto_W3110_Mu_SGS.BroadPeak",
                  'IHFB': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\IHFB_sites_mid_exp_Prieto_W3110_Mu_SGS.BroadPeak",
                  'IHFAB': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\IHFAB_sites_mid_exp_Prieto_W3110_Mu_SGS.BroadPeak",
                  'H-NS': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\H-NS_sites_ME_Kahramanoglou_W3110_Mu_SGS.BroadPeak",
                  'Fis': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Fis_sites_ME_Kahramanoglou_W3110_Mu_SGS.BroadPeak",
                  'MatP': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MatP_sites_37deg_Nolivos_W3110_Mu_SGS.BroadPeak",
                  'MatS': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MatS_sites_Mercier_W3110_Mu_SGS.BroadPeak",
                  'TopoIV' : "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\TopoIV_trusted_sites_Sayyed_W3110_Mu_SGS.BroadPeak",
                  'MukB' : "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MukB_sites_37deg_Nolivos_W3110_Mu_SGS.BroadPeak",
                  'Macrodomains': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Macrodomains_E_coli_w3110_edited.broadPeak"}

path_to_intervals_sets={'MatP 0': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MatP_sites_37deg_Nolivos_W3110_Mu_SGS.BroadPeak",
                        'MatP 3': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MatP_sites_thr_3_37deg_Nolivos_W3110_Mu_SGS.BroadPeak",
                        'MukB 0' : "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MukB_sites_37deg_Nolivos_W3110_Mu_SGS.BroadPeak",
                        'MukB 2' : "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MukB_sites_thr_2_37deg_Nolivos_W3110_Mu_SGS.BroadPeak",
                        'MukB 3' : "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MukB_sites_thr_3_37deg_Nolivos_W3110_Mu_SGS.BroadPeak"}
'''
path_to_intervals_sets={'MatP Ter': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MatP_sites_thr_3_37deg_Nolivos_W3110_Mu_SGS_Ter.BroadPeak",
                        'MatP not Ter': "C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\MatP_sites_thr_3_37deg_Nolivos_W3110_Mu_SGS_no_Ter.BroadPeak"}
#Path for intervals analysis output.
Intervals_analysis_outpath="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_association_with_MatP_MukB_dif_thr\\"
if not os.path.exists(Intervals_analysis_outpath):
    os.makedirs(Intervals_analysis_outpath)
 
 
'''Test for ratio of Poisson intensities in two independent samples

Author: Josef Perktold
License: BSD-3

''' 

#######
#Start of the Poisson test
 
#Copied from statsmodels.stats.weightstats
def _zstat_generic2(value, std_diff, alternative):
    '''generic (normal) z-test to save typing

    can be used as ztest based on summary statistics
    '''
    zstat = value / std_diff
    if alternative in ['two-sided', '2-sided', '2s']:
        pvalue = stats.norm.sf(np.abs(zstat))*2
    elif alternative in ['larger', 'l']:
        pvalue = stats.norm.sf(zstat)
    elif alternative in ['smaller', 's']:
        pvalue = stats.norm.cdf(zstat)
    else:
        raise ValueError('invalid alternative')
    return zstat, pvalue


def poisson_twosample(count1, exposure1, count2, exposure2, ratio_null=1,
                      method='score', alternative='2-sided'):
    '''test for ratio of two sample Poisson intensities

    If the two Poisson rates are g1 and g2, then the Null hypothesis is

    H0: g1 / g2 = ratio_null

    against one of the following alternatives

    H1_2-sided: g1 / g2 != ratio_null
    H1_larger: g1 / g2 > ratio_null
    H1_smaller: g1 / g2 < ratio_null

    Parameters
    ----------
    count1: int
        Number of events in first sample
    exposure1: float
        Total exposure (time * subjects) in first sample
    count2: int
        Number of events in first sample
    exposure2: float
        Total exposure (time * subjects) in first sample
    ratio: float
        ratio of the two Poisson rates under the Null hypothesis. Default is 1.
    method: string
        Method for the test statistic and the p-value. Defaults to `'score'`.
        Current Methods are based on Gu et. al 2008
        Implemented are 'wald', 'score' and 'sqrt' based asymptotic normal
        distribution, and the exact conditional test 'exact-cond', and its mid-point
        version 'cond-midp', see Notes
    alternative : string
        The alternative hypothesis, H1, has to be one of the following

           'two-sided': H1: ratio of rates is not equal to ratio_null (default)
           'larger' :   H1: ratio of rates is larger than ratio_null
           'smaller' :  H1: ratio of rates is smaller than ratio_null

    Returns
    -------
    stat, pvalue two-sided

    not yet
    #results : Results instance
    #    The resulting test statistics and p-values are available as attributes.


    Notes
    -----
    'wald': method W1A, wald test, variance based on separate estimates
    'score': method W2A, score test, variance based on estimate under Null
    'wald-log': W3A
    'score-log' W4A
    'sqrt': W5A, based on variance stabilizing square root transformation
    'exact-cond': exact conditional test based on binomial distribution
    'cond-midp': midpoint-pvalue of exact conditional test

    The latter two are only verified for one-sided example.

    References
    ----------
    Gu, Ng, Tang, Schucany 2008: Testing the Ratio of Two Poisson Rates,
    Biometrical Journal 50 (2008) 2, 2008

    '''
    
    if count1==0 or count2==0:
        return ('NA', 'NA')

    # shortcut names
    y1, n1, y2, n2 = count1, exposure1, count2, exposure2
    d = n2 / n1
    r = ratio_null
    r_d = r / d

    if method in ['score']:
        stat = (y1 - y2 * r_d) / np.sqrt((y1 + y2) * r_d)
        dist = 'normal'
    elif method in ['wald']:
        stat = (y1 - y2 * r_d) / np.sqrt(y1 + y2 * r_d**2)
        dist = 'normal'
    elif method in ['sqrt']:
        stat = 2 * (np.sqrt(y1 + 3 / 8.) - np.sqrt((y2 + 3 / 8.) * r_d))
        stat /= np.sqrt(1 + r_d)
        dist = 'normal'
    elif method in ['exact-cond', 'cond-midp']:
        from statsmodels.stats import proportion
        bp = r_d / (1 + r_d)
        y_total = y1 + y2
        stat = None
        pvalue = proportion.binom_test(y1, y_total, prop=bp, alternative=alternative)
        if method in ['cond-midp']:
            # not inplace in case we still want binom pvalue
            pvalue = pvalue - 0.5 * stats.binom.pmf(y1, y_total, bp)

        dist = 'binomial'

    if dist == 'normal':
        return _zstat_generic2(stat, 1, alternative)
    else:
        return stat, pvalue    

#End of the Poisson test
#######    

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
                tu['Transcription level']=float(line[5].replace(',','.')) #Transcription level
                tu['Description']=str(line[6]) #Description
                ar.append(tu)
                if line[4]=='+':
                    plus+=1
                elif line[4]=="-":
                    minus+=1
        filein.close()
        TUs_sets[k]=ar
        print('Number of TUs in forward for ' + str(k) + ' set: ' + str(plus))
        print('Number of TUs in reverse for ' + str(k) + ' set: ' + str(minus))
    return TUs_sets

#######
#GCSs association with TUs.
#######

def TU_association(GCSs_sets_dict, TU_set, set_type, window_width, path_out, set_name):
    fileout=open(path_out + set_name + '_numbers_of_associated_GCSs.txt', 'w')
    if set_type=='16S operons':
        fileout.write("Operon_ID\tGenes_ID\tStart\tEnd\tStrand\tExpression\t") 
        for k in range(len(GCSs_sets_dict)):
            fileout.write("Condition\tGCSs in US\tGCSs in GB\tGCSs in DS\t")
        fileout.write("\n")
    elif set_type=='operons':
        fileout.write("Operon_ID\tGenes_ID\tStart\tEnd\tStrand\tExpression\t") 
        for k in range(len(GCSs_sets_dict)):
            fileout.write("Condition\tGCSs in USUS\tGCSs in USGB\tGCSs in GBDS\tGCSs in DSDS\t")
        fileout.write("\n")        
    elif set_type=='genes':
        fileout.write("Gene_ID\tGene_name\tStart\tEnd\tStrand\tExpression\t")  
        for k in range(len(GCSs_sets_dict)):
            fileout.write("Condition\tGCSs in USUS\tGCSs in USGB\tGCSs in GBDS\tGCSs in DSDS\t")
        fileout.write("\n")
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
                    operons_stat_ds[a]=[0, [], [], 0, [], [], 0, [], []] #[number of GCSs, N3E values, score values] for US, GB, DS correspondingly.
                elif set_type=='operons' or set_type=='genes':
                    operons_stat_ds[a]=[0, [], [], 0, [], [], 0, [], [], 0, [], []] #[number of GCSs, N3E values, score values] for USUS, USGB, GBDS, DSDS correspondingly.
            fileout.write(a + '\t')
            stats=[0, 0, 0, 0, 0] #USUS USGB GB GBDS DSDS
            N3E=[[], [], [], [], []] #USUS USGB GB GBDS DSDS
            score=[[], [], [], [], []] #USUS USGB GB GBDS DSDS
            for k, v in s.items(): #k - GCS coordinate, v - N3E and score.
                if j['Start']+window_width>k>j['Start'] and j['Strand']=='-': #GBDS
                    stats[3]+=1
                    N3E[3].append(v[0])
                    score[3].append(v[1])                    
                if j['Start']+window_width>k>j['Start'] and j['Strand']=='+': #USGB
                    stats[1]+=1
                    N3E[1].append(v[0])
                    score[1].append(v[1])
                if j['End']>k>j['End']-window_width and j['Strand']=='-': #USGB
                    stats[1]+=1
                    N3E[1].append(v[0])
                    score[1].append(v[1])  
                if j['End']>k>j['End']-window_width and j['Strand']=='+': #GBDS
                    stats[3]+=1
                    N3E[3].append(v[0])
                    score[3].append(v[1])                 
                if j['Start']>k>j['Start']-window_width and j['Strand']=='-': #DSDS
                    stats[4]+=1
                    N3E[4].append(v[0])
                    score[4].append(v[1])
                if j['Start']>k>j['Start']-window_width and j['Strand']=='+': #USUS
                    stats[0]+=1   
                    N3E[0].append(v[0])
                    score[0].append(v[1])
                if j['End']+window_width>k>j['End'] and j['Strand']=='-': #USUS
                    stats[0]+=1
                    N3E[0].append(v[0])
                    score[0].append(v[1])
                if j['End']+window_width>k>j['End'] and j['Strand']=='+': #DSDS
                    stats[4]+=1 
                    N3E[4].append(v[0])
                    score[4].append(v[1])
                if j['End']>k>j['Start']: #GB
                    stats[2]+=1
                    N3E[2].append(v[0])
                    score[2].append(v[1])
            if set_type=='16S operons':
                fileout.write(str(stats[0]) + '\t' + str(stats[2]) + '\t' + str(stats[4]) + '\t') #US, GB, DS
                #US
                operons_stat_ds[a][0]+=stats[0]
                operons_stat_ds[a][1]+=N3E[0]
                operons_stat_ds[a][2]+=score[0] 
                #GB
                operons_stat_ds[a][3]+=stats[2]
                operons_stat_ds[a][4]+=N3E[2]
                operons_stat_ds[a][5]+=score[2] 
                #DS
                operons_stat_ds[a][6]+=stats[4]
                operons_stat_ds[a][7]+=N3E[4]
                operons_stat_ds[a][8]+=score[4]                 
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

def TU_interval_stat_analysis(GCSs_sets_dict, intervals_GCSs_dict, intervals, score_data, window_width, set_type, path_out, set_name, genome_len):
    #Correcting genome length.
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    del_len=0
    for i in deletions:
        del_len+=i[1]-i[0]
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
            
    fileout=open(path_out + set_name + '_compartments_statistics_GCSs_number_N3E_score.txt', 'w')
    #Performes t-test for comparison of intervals GCSs N3E and scores with all GCSs N3E and scores.
    #Performes binomail test for enrichment estimation of GCSs fall into intervals.
    fileout.write('Test\tAntibiotic\t')
    for i in range(len(intervals_GCSs_dict.values())-1):
        fileout.write('Examined value\tValue (intervals)\tValue (overall)\tp-value\tAdditional\t')
    fileout.write('\n') 
    
    if set_type=='16S operons':
        for a, ns in GCSs_values_dict.items():
            #N3E
            N3E_stat=poisson_twosample(len(ns[0]), sum(ns[0]), len(intervals_GCSs_dict[a][1]), sum(intervals_GCSs_dict[a][1]), ratio_null=1,
                                  method='exact-cond', alternative='two-sided') #N3E US           
            fileout.write('t-test\t' + a + '\tUS N3E\t' + str(round(np.mean(intervals_GCSs_dict[a][1]),3)) + '\t' + 
                          str(round(np.mean(ns[0]),3)) + '\t' + str(N3E_stat[1]) + '\n') 
            N3E_stat=poisson_twosample(len(ns[0]), sum(ns[0]), len(intervals_GCSs_dict[a][4]), sum(intervals_GCSs_dict[a][4]), ratio_null=1,
                                  method='exact-cond', alternative='two-sided') #N3E GB             
            fileout.write('t-test\t' + a + '\tGB N3E\t' + str(round(np.mean(intervals_GCSs_dict[a][4]),3)) + '\t' + 
                          str(round(np.mean(ns[0]),3)) + '\t' + str(N3E_stat[1]) + '\n') 
            N3E_stat=poisson_twosample(len(ns[0]), sum(ns[0]), len(intervals_GCSs_dict[a][7]), sum(intervals_GCSs_dict[a][7]), ratio_null=1,
                                  method='exact-cond', alternative='two-sided') #N3E DS             
            fileout.write('t-test\t' + a + '\tDS N3E\t' + str(round(np.mean(intervals_GCSs_dict[a][7]),3)) + '\t' + 
                          str(round(np.mean(ns[0]),3)) + '\t' + str(N3E_stat[1]) + '\n') 
            #Score
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][2]) #Score US
            fileout.write('t-test\t' + a + '\tUS Score\t' + str(round(np.mean(intervals_GCSs_dict[a][2]),3)) + '\t' + 
                          str(round(np.mean(ns[1]),3)) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(round(score_stat[0],3)) + '\n')
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][5]) #Score GB
            fileout.write('t-test\t' + a + '\tGB Score\t' + str(round(np.mean(intervals_GCSs_dict[a][5]),3)) + '\t' + 
                          str(round(np.mean(ns[1]),3)) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(round(score_stat[0],3)) + '\n')
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][8]) #Score DS
            fileout.write('t-test\t' + a + '\tDS Score\t' + str(round(np.mean(intervals_GCSs_dict[a][8]),3)) + '\t' + 
                          str(round(np.mean(ns[1]),3)) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(round(score_stat[0],3)) + '\n')
            #Number of GCSs
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][0], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in US
            fileout.write('binomial test\t' + a + '\tUS Number of GCSs\t' + str(intervals_GCSs_dict[a][0]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\n')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][3], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in GB
            fileout.write('binomial test\t' + a + '\tGB Number of GCSs\t' + str(intervals_GCSs_dict[a][3]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\n')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][6], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in DS
            fileout.write('binomial test\t' + a + '\tDS Number of GCSs\t' + str(intervals_GCSs_dict[a][6]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\n')            
    
    elif set_type=='operons' or set_type=='genes':
        for a, ns in GCSs_values_dict.items():
            #N3E
            N3E_stat=poisson_twosample(len(ns[0]), sum(ns[0]), len(intervals_GCSs_dict[a][1]), sum(intervals_GCSs_dict[a][1]), ratio_null=1,
                                  method='exact-cond', alternative='two-sided') #N3E USUS              
            fileout.write('t-test\t' + a + '\tUSUS N3E\t' + str(round(np.mean(intervals_GCSs_dict[a][1]),3)) + '\t' + 
                          str(round(np.mean(ns[0]),3)) + '\t' + str(N3E_stat[1]) + '\t')
            N3E_stat=poisson_twosample(len(ns[0]), sum(ns[0]), len(intervals_GCSs_dict[a][4]), sum(intervals_GCSs_dict[a][4]), ratio_null=1,
                                  method='exact-cond', alternative='two-sided') #N3E USGB              
            fileout.write('USGB N3E\t' + str(round(np.mean(intervals_GCSs_dict[a][4]),3)) + '\t' + 
                          str(round(np.mean(ns[0]),3)) + '\t' + str(N3E_stat[1]) + '\t')               
            N3E_stat=poisson_twosample(len(ns[0]), sum(ns[0]), len(intervals_GCSs_dict[a][7]), sum(intervals_GCSs_dict[a][7]), ratio_null=1,
                                  method='exact-cond', alternative='two-sided') #N3E GBDS
            fileout.write('GBDS N3E\t' + str(round(np.mean(intervals_GCSs_dict[a][7]),3)) + '\t' + 
                          str(round(np.mean(ns[0]),3)) + '\t' + str(N3E_stat[1]) + '\t')   
            N3E_stat=poisson_twosample(len(ns[0]), sum(ns[0]), len(intervals_GCSs_dict[a][10]), sum(intervals_GCSs_dict[a][10]), ratio_null=1,
                                  method='exact-cond', alternative='two-sided') #N3E DSDS
            fileout.write('DSDS N3E\t' + str(round(np.mean(intervals_GCSs_dict[a][10]),3)) + '\t' + 
                          str(round(np.mean(ns[0]),3)) + '\t' + str(N3E_stat[1]) + '\n')               
            #Score
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][2]) #Score USUS
            fileout.write('t-test\t' + a + '\tUSUS Score\t' + str(round(np.mean(intervals_GCSs_dict[a][2]),3)) + '\t' + 
                          str(round(np.mean(ns[1]),3)) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(round(score_stat[0],3)) + '\t')
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][5]) #Score USGB
            fileout.write('USGB Score\t' + str(round(np.mean(intervals_GCSs_dict[a][5]),3)) + '\t' + 
                          str(round(np.mean(ns[1]),3)) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(round(score_stat[0],3)) + '\t')
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][8]) #Score GBDS
            fileout.write('GBDS Score\t' + str(round(np.mean(intervals_GCSs_dict[a][8]),3)) + '\t' + 
                          str(round(np.mean(ns[1]),3)) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(round(score_stat[0],3)) + '\t')
            score_stat=stats.ttest_ind(ns[1], intervals_GCSs_dict[a][11]) #Score DSDS
            fileout.write('DSDS Score\t' + str(round(np.mean(intervals_GCSs_dict[a][11]),3)) + '\t' + 
                          str(round(np.mean(ns[1]),3)) + '\t' + str(score_stat[1]) + '\t' + 't-statistic: ' + str(round(score_stat[0],3)) + '\n')            
            #Number of GCSs
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][0], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in USUS
            fileout.write('binomial test\t' + a + '\tUSUS Number of GCSs\t' + str(intervals_GCSs_dict[a][0]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\tNA\t')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][3], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in USGB
            fileout.write('USGB Number of GCSs\t' + str(intervals_GCSs_dict[a][3]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\tNA\t')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][6], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in GBDS
            fileout.write('GBDS Number of GCSs\t' + str(intervals_GCSs_dict[a][6]) + '\t' + 
                          str(len(ns[0])) + '\t' + str(GCSs_number_stat) + '\tNA\t')
            GCSs_number_stat=binom.cdf(intervals_GCSs_dict[a][9], len(ns[0]), intervals_GCSs_dict[a][-1]/genome_len_dc) #Number of GCSs in DSDS
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
        fileout.write('t-test\t' + 'NA' + '\t DS Interval score\t' + str(round(np.mean(intervals_values),3)) + '\t' + 
                      str(round(np.mean(score_data),3)) + '\t' + str(intervals_score_stat[1]) + '\t' + 't-statistic: ' + 
                      str(round(intervals_score_stat[0],3)) + '\n')
    
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
                      '\t USUS Interval score\t' + str(round(np.mean(intervals_values[0]),3)) + '\t' + str(round(np.mean(score_data),3)) + '\t' + 
                      str(intervals_score_stat[0][1]) + '\t' + 't-statistic: ' + str(round(intervals_score_stat[0][0],3)) + '\t' + 
                      '\t USGB Interval score\t' + str(round(np.mean(intervals_values[1]),3)) + '\t' + str(round(np.mean(score_data),3)) + '\t' + 
                      str(intervals_score_stat[1][1]) + '\t' + 't-statistic: ' + str(round(intervals_score_stat[1][0],3)) + '\t' + 
                      '\t GBDS Interval score\t' + str(round(np.mean(intervals_values[2]),3)) + '\t' + str(round(np.mean(score_data),3)) + '\t' + 
                      str(intervals_score_stat[2][1]) + '\t' + 't-statistic: ' + str(round(intervals_score_stat[2][0],3)) + '\t' + 
                      '\t DSDS Interval score\t' + str(round(np.mean(intervals_values[3]),3)) + '\t' + str(round(np.mean(score_data),3)) + '\t' + 
                      str(intervals_score_stat[3][1]) + '\t' + 't-statistic: ' + str(round(intervals_score_stat[3][0],3)) + '\n')      
    fileout.close() 
    return

#######
#Statistical analysis for the number of GCSs associated with TUs, additional normalization step is added. 
#######

def GCSs_num_normalization(all_TUs, part_TUs, GCSs_sets_dict, all_TUs_set, part_TUs_set):
    GCSs_set_exp_interval_dict={} 
    for a, s in all_TUs.items():
        GCSs_all=all_TUs[a][0]+all_TUs[a][3]+all_TUs[a][6]+all_TUs[a][9] #Total number of GCSs observed to be associated with the full set of TUs.
        GCSs_all_exp_compartment=float(GCSs_all)/4 #Expected number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for full set of TUs.
        GCSs_set_exp_interval=GCSs_all_exp_compartment*len(part_TUs_set)/len(all_TUs_set) #Expected number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for partial set of TUs.
        GCSs_set_exp_interval_dict[a]=[GCSs_set_exp_interval] 
        for j in range(int((len(s)-1)/3)):
            GCSs_set_exp_interval_dict[a].append(part_TUs[a][j*3]) #Number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for partial set of TUs.
            GCSs_set_exp_interval_dict[a].append(binom.cdf(part_TUs[a][j*3], GCSs_set_exp_interval*4, 0.25)) #p-value of binomial test for the number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for partial set of TUs.
            GCSs_set_exp_interval_dict[a].append(float(part_TUs[a][j*3])*100000/(len(GCSs_sets_dict[a])*len(part_TUs_set))) #Normalized number of GCSs fall into particular compartment (USUS, USGB, GBDS, DSDS) for partial set of TUs.
    return GCSs_set_exp_interval_dict #[GCSs expected] + [GCSs obs, p-value, GCSs norm]*[USUS, USGB, GBDS, DSDS]


#######
#Statistical analysis for the number of GCSs associated with 16S operons, additional normalization step is added. 
#######

def GCSs_number_norm_16S(intervals_GCSs_dict, GCSs_sets_dict, genome_len):
    #Correcting genome length.
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    del_len=0
    for i in deletions:
        del_len+=i[1]-i[0]
    genome_len_dc=genome_len-del_len  
    #Number of GCSs normalization and statistics.
    GCSs_set_exp_interval_dict={}
    for a, ns in intervals_GCSs_dict.items():
        Num_GCSs_expected=len(GCSs_sets_dict[a])*ns[-1]/genome_len_dc
        GCSs_set_exp_interval_dict[a]=[Num_GCSs_expected] 
        for j in range(int((len(ns)-1)/3)):
            GCSs_set_exp_interval_dict[a].append(ns[j*3]) #Number of GCSs fall into particular compartment (US, GB, DS) of 16S operons.
            GCSs_set_exp_interval_dict[a].append(binom.cdf(ns[j*3], len(GCSs_sets_dict[a]), ns[-1]/genome_len_dc)) #p-value of binomial test for the number of GCSs fall into particular compartment (US, GB, DS) of the 16S operons.
            GCSs_set_exp_interval_dict[a].append(float(ns[j*3])/Num_GCSs_expected) #Normalized number of GCSs fall into particular compartment (US, GB, DS) of the 16S operons.
    return GCSs_set_exp_interval_dict #[GCSs expected] + [GCSs obs, p-value, GCSs norm]*[US, GB, DS]


#######
#Writes GCSs numbers information to file: Number of GCSs expected, Number of GCSs observed, p-value - binomial test, Number of GCSs normalized.
#######

def write_GCSs_norm(GCSs_set_exp_interval_dict, path_out, set_name):
    fileout=open(path_out + set_name + '_normalized_GCSs_numbers_and_statistics.txt', 'w')
    fileout.write('Condition\tCompartment\tNumber of GCSs expected\tNumber of GCSs observed\tp-value\tNumber of GCSs normalized\n')
    if set_name=='16S_operons':
        Compartment_names=['US', 'GB', 'DS']
    else:
        Compartment_names=['USUS', 'USGB', 'GBDS', 'DSDS']
    for a, s in GCSs_set_exp_interval_dict.items():
        for i in range(len(Compartment_names)):
            fileout.write(a + '\t' + Compartment_names[i] + '\t' + str(round(s[0],3)) + '\t' + str(s[(i*3)+1]) + '\t' + str(s[(i*3)+2]) + '\t' + str(round(s[(i*3)+3],3)) +'\n')
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
        if k in ['Macrodomains']:
            for line in filein:
                line=line.rstrip().split('\t')
                int_start=int(line[1])
                int_end=int(line[2])
                ar.append([int_start, int_end, line[3]]) #Macrodomain start, macrodomain end, macrodomain name
        else:
            for line in filein:
                line=line.rstrip().split('\t')
                int_start=int(line[1])
                int_end=int(line[2])
                del_check=0
                for j in range(len(deletions)):
                    if deletions[j][1]>int_start>deletions[j][0] or deletions[j][1]>int_end>deletions[j][0]:
                        del_check=1
                if del_check==0:
                    ar.append([int_start, int_end]) #Interval start, interval end        
        intervals_sets_dict[k]=ar
        print("Number of " + str(k) + " regions: " + str(len(ar)))
        filein.close()
    return intervals_sets_dict

#######
#GCSs in intervals.
#######   

def GCSs_in_intervals(GCSs_sets_dict, intervals, score_data, path_out, genome_len):
    #Preparing output files.
    fileout=open(path_out+'GCSs_associated_with_intervals_statistics.txt', 'w') #For all intervals sets.
    fileout.write('Interval type\tCondition\tNumber of GCSs expected\tNumber of GCSs observed\tp-value\t' +
                  'GCSs N3E (intervals)\tGCSs N3E (overall)\tp-value\t' + 
                  'GCSs score (intervals)\tGCSs score (overall)\tp-value\tt-statistic\n')
    
    fileout1=open(path_out+'GCSs_association_with_BIME.txt', 'w') #Specially for BIMEs1, BIMEs2 sets.
    fileout1.write('Interval type\tStart\tEnd\tCfx (number of GCSs)\tRifCfx (number of GCSs)\tMicro (number of GCSs)\tOxo (number of GCSs)\n')   
    fileout_bime_rel=open(path_out+'GCSs_relative_positions_within_BIME.txt', 'w')
    fileout_bime_rel.write('Interval type\tStart\tEnd\tInterval len\tCfx GCSs relative coordinates\tRifCfx GCSs relative coordinates\tMicro GCSs relative coordinates\tOxo GCSs relative coordinates\n')
    fileout_macro=open(path_out+'GCSs_association_with_macrodomains.txt', 'w') #Specially for chromosomal macrodomains.
    fileout_macro.write('Macrodomain\tStart\tEnd')   
    for i in range(len(GCSs_sets_dict)):
        fileout_macro.write('\tCondition\tNumber of GCSs observed\tNumber of GCSs expected\tp-value (binomial test)')
    fileout_macro.write('\n')
    
    #Correcting genome length.
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]] #Deletions in E. coli w3110 strain that corresponds to DY330 strain.
    del_len=0
    for i in deletions:
        del_len+=i[1]-i[0]
    genome_len_dc=genome_len-del_len  
    
    #Prepares lists of GCSs N3E and scores for statistical tests.
    GCSs_values_dict={}
    for a, s in GCSs_sets_dict.items():
        for k, v in s.items():
            if a not in GCSs_values_dict:
                GCSs_values_dict[a]=[[v[0]], [v[1]]]
            else:
                GCSs_values_dict[a][0].append(v[0]) #N3E 
                GCSs_values_dict[a][1].append(v[1]) #Score
    
    #Finds GCSs associated with intervals.
    GCSs_associated_info={} #Dictionary contains element corresponds to interval sets.
    for k, v in intervals.items(): #Iterate different interval types.
        GCSs_associated_info[k]={}
        intervals_len=0
        for i in v: #Iterate intevals.
            intervals_len+=i[1]-i[0]  
            interval_associated_GCSs={}
            for a, s in GCSs_sets_dict.items(): #Iterate GCSs sets.
                if a not in GCSs_associated_info[k]:
                    GCSs_associated_info[k][a]=[0, [], []] #Dictionary contains elements corresponds to GCSs sets. [Number of GCSs, N3E values, Score values]
                interval_associated_GCSs[a]=[] #Counter for the number of GCSs fall into the particular interval. Stores GCSs coordinates.
                for gcs, info in s.items(): #Iterate GCSs.
                    if i[1]>=i[0]>=0: #Interval does not cross position for genome start/end
                        if i[1]>gcs>i[0]:
                            GCSs_associated_info[k][a][0]+=1
                            GCSs_associated_info[k][a][1].append(info[0])
                            GCSs_associated_info[k][a][2].append(info[1])
                            interval_associated_GCSs[a].append(gcs)
                    elif 0<i[1]<i[0]: #Interval crosses position for genome start/end
                        if i[1]>gcs>=0 or genome_len>gcs>=i[0]:
                            GCSs_associated_info[k][a][0]+=1
                            GCSs_associated_info[k][a][1].append(info[0])
                            GCSs_associated_info[k][a][2].append(info[1])
                            interval_associated_GCSs[a].append(gcs)                             
                    elif i[0]<0 and i[1]>0: #Interval crosses position for genome start/end (alternative coordinate system for interval boundary setting)
                        if i[1]>gcs>=0 or genome_len>gcs>=genome_len+i[0]:
                            GCSs_associated_info[k][a][0]+=1
                            GCSs_associated_info[k][a][1].append(info[0])
                            GCSs_associated_info[k][a][2].append(info[1])
                            interval_associated_GCSs[a].append(gcs)                          
            
            #Writing data for BIMEs
            if k in ['BIMEs1', 'BIMEs2']:
                fileout1.write(k + '\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + str(len(interval_associated_GCSs['Cfx'])) + '\t' + 
                               str(len(interval_associated_GCSs['RifCfx'])) + '\t' + str(len(interval_associated_GCSs['Micro'])) + '\t' +
                               str(len(interval_associated_GCSs['Oxo'])) + '\n')  
                #Check whether no GCSs fall into this particular interval.
                Num_in_int=0
                for cond, gcs_nums in interval_associated_GCSs.items():
                    Num_in_int+=len(gcs_nums)
                #If GCSs were found, then write.
                if Num_in_int!=0:
                    fileout_bime_rel.write(k + '\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + str(i[1]-i[0]))
                    #Cfx GCSs
                    if len(interval_associated_GCSs['Cfx'])==0:
                        fileout_bime_rel.write('\tNA')
                    else:
                        coords='\t'
                        for gcs in interval_associated_GCSs['Cfx']:
                            coords+=str(gcs-i[0])+':'
                        fileout_bime_rel.write(coords.rstrip(':'))
                    #RifCfx GCSs
                    if len(interval_associated_GCSs['RifCfx'])==0:
                        fileout_bime_rel.write('\tNA')
                    else:
                        coords='\t'
                        for gcs in interval_associated_GCSs['RifCfx']:
                            coords+=str(gcs-i[0])+':'
                        fileout_bime_rel.write(coords.rstrip(':'))
                    #Micro GCSs
                    if len(interval_associated_GCSs['Micro'])==0:
                        fileout_bime_rel.write('\tNA')
                    else:
                        coords='\t'
                        for gcs in interval_associated_GCSs['Micro']:
                            coords+=str(gcs-i[0])+':'
                        fileout_bime_rel.write(coords.rstrip(':'))
                    #Oxo GCSs
                    if len(interval_associated_GCSs['Oxo'])==0:
                        fileout_bime_rel.write('\tNA')
                    else:
                        coords='\t'
                        for gcs in interval_associated_GCSs['Oxo']:
                            coords+=str(gcs-i[0])+':'
                        fileout_bime_rel.write(coords.rstrip(':'))
                    fileout_bime_rel.write('\n')
                
            #Writing data for chromosomal macrodomains       
            elif k in ['Macrodomains']:
                fileout_macro.write(i[2] + '\t' + str(i[0]) + '\t' + str(i[1]))
                macro_len=i[1]-i[0]
                for delet in deletions:
                    where_is_del=0
                    if i[1]>delet[0]>i[0]: #Start of the deletion is in the macrodomain
                        where_is_del=1
                    if i[1]>delet[1]>i[0]: #End of the deletion is in the macrodomain
                        where_is_del+=2
                    if delet[1]>i[1]>i[0]>delet[0]: #Macrodomain is deleted
                        macro_len=0
                        continue
                    if where_is_del==0: #Hence particular deletion does not interrupt the macrodomain
                        continue
                    elif where_is_del==1: #Hence deletion covers the end of the macrodomain
                        macro_len+=delet[0]-i[1]
                    elif where_is_del==2: #Hence deletion covers the start of the macrodomain
                        macro_len+=i[0]-delet[1]
                    elif where_is_del==3: #Hence macrodomain covers the deletion
                        macro_len+=delet[0]-delet[1]
                for ab, gcs_num in interval_associated_GCSs.items():
                    condition=ab
                    num_GCSs_observed=len(interval_associated_GCSs[ab])
                    num_GCSs_expected=len(GCSs_sets_dict[ab])*macro_len/genome_len_dc
                    num_GCSs_p_value=binom.cdf(num_GCSs_observed, len(GCSs_sets_dict[ab]), macro_len/genome_len_dc)
                    fileout_macro.write('\t' + condition + '\t' + str(num_GCSs_observed) + '\t' + str(round(num_GCSs_expected, 3)) + '\t' + str(num_GCSs_p_value))  
                fileout_macro.write('\n')
        
        #Writing statistics for all kinds of interval except chromosomal macrodomains.
        if k not in ['Macrodomains']:
            for a, s in GCSs_associated_info[k].items():
                relative_intervals_len=float(intervals_len)/genome_len_dc
                total_number_of_GCSs=len(GCSs_sets_dict[a])
                expected_number_of_GCSs_in_interval=relative_intervals_len*total_number_of_GCSs 
                GCSs_num_stat=binom.cdf(s[0], total_number_of_GCSs, relative_intervals_len) #Number of GCSs stat.
                GCSs_N3E_stat=poisson_twosample(len(s[1]), sum(s[1]), len(GCSs_values_dict[a][0]), sum(GCSs_values_dict[a][0]), ratio_null=1,
                                      method='exact-cond', alternative='two-sided') #GCSs N3E stat.                
                GCSs_score_stat=stats.ttest_ind(s[2], GCSs_values_dict[a][1]) #GCSs score stat.
                fileout.write(k + '\t' + a + '\t' + str(round(expected_number_of_GCSs_in_interval,3)) + '\t' + str(s[0]) + '\t' + str(GCSs_num_stat) + '\t' + 
                              str(round(np.mean(s[1]),3)) + '\t' + str(round(np.mean(GCSs_values_dict[a][0]),3)) + '\t' + str(GCSs_N3E_stat[1]) + '\t' + 
                              str(round(np.mean(s[2]),3)) + '\t' + str(round(np.mean(GCSs_values_dict[a][1]),3)) + '\t' + str(GCSs_score_stat[1]) + '\t' + str(round(GCSs_score_stat[0],3)) + '\n')
    fileout.close()
    fileout1.close()
    fileout_macro.close()
    fileout_bime_rel.close()
    
    #Intervals score statistics.
    fileout2=open(path_out+'Intervals_score_statistics.txt', 'w')
    fileout2.write('Test\tInterval\tScore (intervals)\tScore (overall)\tp-value\tt-statistic\n')
    for k, v in intervals.items():
        if k not in ['Macrodomains']:
            intervals_values=[]
            for i in v:
                intervals_values+=score_data[i[0]:i[1]]
            intervals_score_stat=stats.ttest_ind(intervals_values, score_data) #Score
            fileout2.write('t-test\t' + k + '\t' + str(round(np.mean(intervals_values),3)) + '\t' + 
                      str(round(np.mean(score_data),3)) + '\t' + str(intervals_score_stat[1]) + '\t' + 't-statistic: ' + 
                      str(round(intervals_score_stat[0],3)) + '\n')  
    fileout2.close()
    return GCSs_associated_info


#######
#Wrapper for TUs analysis functions.
####### 

def TU_analysis_wrapper(input_dict, inpath, TUs_sets_path, path_out, genome_len):
    #Reading input.
    GCSs_sets_dict=trusted_GCSs_parsing(input_dict) #Parsing GCSs
    score_data=score_data_parser(inpath, 'score') #Parsing score file
    TU_sets_dict=TUs_parser(TUs_sets_path) #Parsing TUs
    
    #16S operons analysis.
    window_width_16S_operons=5000
    set_type_16S='16S operons'
    GCSs_16S_assoc_info=TU_association(GCSs_sets_dict, TU_sets_dict['16S_operons'], set_type_16S, window_width_16S_operons, path_out, set_type_16S)
    TU_interval_stat_analysis(GCSs_sets_dict, GCSs_16S_assoc_info, TU_sets_dict['16S_operons'], score_data, window_width_16S_operons, set_type_16S, path_out, set_type_16S, genome_len)
    GCSs_set_exp_interval_dict_16S=GCSs_number_norm_16S(GCSs_16S_assoc_info, GCSs_sets_dict, genome_len)
    write_GCSs_norm(GCSs_set_exp_interval_dict_16S, path_out, '16S_operons')  
    
    #All genes analysis.
    window_width=650
    GCSs_all_genes_assoc_info=TU_association(GCSs_sets_dict, TU_sets_dict['All_genes'], 'genes', window_width, path_out, 'All_genes')
    TU_interval_stat_analysis(GCSs_sets_dict, GCSs_all_genes_assoc_info, TU_sets_dict['All_genes'], score_data, window_width, 'genes', path_out, 'All_genes', genome_len)
    GCSs_set_exp_interval_dict_ag=GCSs_num_normalization(GCSs_all_genes_assoc_info, GCSs_all_genes_assoc_info, GCSs_sets_dict, TU_sets_dict['All_genes'], TU_sets_dict['All_genes'])
    write_GCSs_norm(GCSs_set_exp_interval_dict_ag, path_out, 'All_genes')    
    
    #Other genes sets analysis.
    for k, v in TU_sets_dict.items():
        if k not in ['All_genes'] and k.find('genes')>0: #k contains 'genes' as a substring but not equial to 'All_genes'.
            GCSs_genes_assoc_info=TU_association(GCSs_sets_dict, v, 'genes', window_width, path_out, k)
            TU_interval_stat_analysis(GCSs_sets_dict, GCSs_genes_assoc_info, v, score_data, window_width, 'genes', path_out, k, genome_len)
            GCSs_set_exp_interval_dict_g=GCSs_num_normalization(GCSs_all_genes_assoc_info, GCSs_genes_assoc_info, GCSs_sets_dict, TU_sets_dict['All_genes'], TU_sets_dict[k])
            write_GCSs_norm(GCSs_set_exp_interval_dict_g, path_out, k)
            
    #All operons analysis.
    window_width=650
    GCSs_all_operons_assoc_info=TU_association(GCSs_sets_dict, TU_sets_dict['All_operons'], 'operons', window_width, path_out, 'All_operons')
    TU_interval_stat_analysis(GCSs_sets_dict, GCSs_all_operons_assoc_info, TU_sets_dict['All_operons'], score_data, window_width, 'operons', path_out, 'All_operons', genome_len)
    GCSs_set_exp_interval_dict_ao=GCSs_num_normalization(GCSs_all_operons_assoc_info, GCSs_all_operons_assoc_info, GCSs_sets_dict, TU_sets_dict['All_operons'], TU_sets_dict['All_operons'])
    write_GCSs_norm(GCSs_set_exp_interval_dict_ao, path_out, 'All_operons')    
    
    #Other operons sets analysis (except 16S operons).
    for k, v in TU_sets_dict.items():
        if k not in ['All_operons', '16S_operons'] and k.find('operons')>0: #k contains 'operons' as a substring but not equial to 'All_operons' or '16S_operons'.
            GCSs_operons_assoc_info=TU_association(GCSs_sets_dict, v, 'operons', window_width, path_out, k)
            TU_interval_stat_analysis(GCSs_sets_dict, GCSs_operons_assoc_info, v, score_data, window_width, 'operons', path_out, k, genome_len)
            GCSs_set_exp_interval_dict_o=GCSs_num_normalization(GCSs_all_operons_assoc_info, GCSs_operons_assoc_info, GCSs_sets_dict, TU_sets_dict['All_operons'], TU_sets_dict[k])
            write_GCSs_norm(GCSs_set_exp_interval_dict_o, path_out, k)
    return


TU_analysis_wrapper(path_to_GCSs_files, Score_path, path_to_TUs_sets, TU_analysis_outpath, Genome_length)

#######
#Wrapper for intervals analysis functions.
####### 

def Interval_analysis_wrapper(input_dict, inpath, intervals_sets_path, path_out, genome_len):
    #Reading input.
    GCSs_sets_dict=trusted_GCSs_parsing(input_dict) #Parsing GCSs
    score_data=score_data_parser(inpath, 'score') #Parsing score file   
    Intervals_sets_dict=broadpeak_pars(intervals_sets_path) #Parsing intervals
    
    #Statistics calculation.
    GCSs_in_intervals(GCSs_sets_dict, Intervals_sets_dict, score_data, path_out, genome_len)
    return

#Interval_analysis_wrapper(path_to_GCSs_files, Score_path, path_to_intervals_sets, Intervals_analysis_outpath, Genome_length)

print('Script ended its work succesfully!') 