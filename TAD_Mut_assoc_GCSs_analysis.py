###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script analysis colocalization of mutations (Foster, 2015) and TADs borders (Lioy, 2018) with GCSs.
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
from scipy.stats import chisquare

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Input data - GCSs, TAB.
path_to_GCSs_files={'Cfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\Cfx_10mkM_trusted_GCSs_h_s.txt",
                    'RifCfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\RifCfx_trusted_GCSs_h_s.txt",
                    'Micro': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\Micro_trusted_GCSs_h_s.txt",
                    'Oxo': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\Oxo_trusted_GCSs_h_s.txt"} 
#Input data - transcription, TAB.
Transcription_data_path='C:\Sutor\science\DNA-gyrase\Results\Final_data_2\Expression_data\Deletion_corrected\DOOR_Mu_del_cor_genes_expression.txt'
#TAD coordinates.
TAD_150_inputpath='C:\Sutor\science\DNA-gyrase\Results\Hi-C\TAD_calling\CIDs_150kb_31_W3110_Mu_compatible.bed'
TAD_80_inputpath='C:\Sutor\science\DNA-gyrase\Results\Hi-C\TAD_calling\CIDs_80kb_57_W3110_Mu_compatible.bed'
#Mutations coordinates.
Mut_wt_inputpath='C:\Sutor\science\E_coli_variome\Foster_data\Foster,2015,WT_E_coli_w3110_Mu_compatible.vcf'
Mut_mut_inputpath='C:\Sutor\science\E_coli_variome\Foster_data\Foster,2015,Mut_E_coli_w3110_Mu_compatible.vcf'
#Output log file
Outputpath="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_association_with_TADs_and_mutations\\"
if not os.path.exists(Outputpath):
    os.makedirs(Outputpath)
Outputpath+="GCSs_association_with_TADs_interTADs_mutations.txt"


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
    print('\n')
    return GCSs_sets_dict

#######
#TAD coordinates parser (BED).
#######

def TAD_parser(TAD_inpath):
    TAD_in=open(TAD_inpath, 'r')
    TAD_coord=[]
    for line in TAD_in:
        line=line.rstrip().split('\t')
        TAD_coord.append([int(line[1])+1, int(line[2])+1]) #[TAD start, TAD end]
    TAD_in.close()
    return TAD_coord

#######
#Mutations parser (VCFv4.1).
#######

def Mutations_parser(Mut_inpath):
    Mut_in=open(Mut_inpath, 'r')
    Mut_dict={}
    for line in Mut_in:
        if line[0] not in ['#']:
            line=line.rstrip().split('\t')
            Mut_dict[int(line[1])]=line
    Mut_in.close()
    return Mut_dict

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
    print('Whole genome average transcription: ' + str(sum(transcription)/len(transcription)))
    return transcription

#######
#Colocalization.
#######

def colocal_search(GCSs_sets_dict, transcription, data_type, data_to_coloc, Bad_TAD_reg, TAD_inpath, Mut_inpath, outpath):
    fileout=open(outpath, 'a')
    genome_len=4647454 #Overall genome len
    win_width=1000
        
    if data_type=='TAD':
        #Genome len correction on bad TAD regions (deletions are included into bad data).
        TAD_genome_len=genome_len
        for i in Bad_TAD_reg:
            TAD_genome_len+=i[0]-i[1] #Genome len corrected on bad TADs
        #Overall len of TADs
        sum_TAD_len=0
        for TAD in data_to_coloc:
            sum_TAD_len+=TAD[1]-TAD[0] #Len of TADs
        #Overall len on interTADs  
        sum_interTAD_len=TAD_genome_len-sum_TAD_len #Len of interTADs
        #GCSs association
        TAD_assoc_dict={}
        for a, gcss_set in GCSs_sets_dict.items(): #Itarates GCSs sets (conditions)
            TAD_assoc_dict[a]={'GCSs in TAD':0, 'GCSs in interTAD':0, 'GCSs in bad TAD':0, 'GCSs in TAD border':0, 'GCSs in interTAD border':0}
            for k, v in gcss_set.items(): #Itarates GCSs
                status=0
                for TAD in data_to_coloc: #Iterates TADs
                    if TAD[1]>=k>=TAD[0]: #GCS falls into TAD
                        TAD_assoc_dict[a]['GCSs in TAD']+=1
                        status=1
                    if TAD[0]+win_width>=k>=TAD[0] or TAD[1]>=k>=TAD[1]-win_width:
                        TAD_assoc_dict[a]['GCSs in TAD border']+=1
                    if TAD[1]+win_width>=k>=TAD[1] or TAD[0]>=k>=TAD[0]-win_width:
                        TAD_assoc_dict[a]['GCSs in interTAD border']+=1
                if status==0: #GCS not in TAD
                    for bad_TAD in Bad_TAD_reg:
                        if bad_TAD[1]>k>bad_TAD[0]: #GCSs fall in the bad TAD region.
                            status=1
                            TAD_assoc_dict[a]['GCSs in bad TAD']+=1
                            break
                    if status==0: #GCSs is in interTAD
                        TAD_assoc_dict[a]['GCSs in interTAD']+=1
            #TADs interTADs defenitions
            GCSs_in_TADs=TAD_assoc_dict[a]['GCSs in TAD']
            GCSs_in_interTADs=TAD_assoc_dict[a]['GCSs in interTAD']
            prob_of_TAD=sum_TAD_len/(sum_interTAD_len+sum_TAD_len)
            prob_of_interTAD=sum_interTAD_len/(sum_interTAD_len+sum_TAD_len)
            GCSs_in_TADs_exp=(GCSs_in_TADs+GCSs_in_interTADs)*prob_of_TAD
            GCSs_in_interTADs_exp=(GCSs_in_TADs+GCSs_in_interTADs)*prob_of_interTAD
            #Borders defenitions
            GCSs_in_TADs_border=TAD_assoc_dict[a]['GCSs in TAD border']
            GCSs_in_interTAD_border=TAD_assoc_dict[a]['GCSs in interTAD border']
            prob_of_TAD_border=2*win_width*len(data_to_coloc)/(sum_interTAD_len+sum_TAD_len)
            prob_of_interTAD_border=prob_of_TAD_border
            GCSs_in_TADs_border_exp=(GCSs_in_TADs+GCSs_in_interTADs)*prob_of_TAD_border
            GCSs_in_interTADs_border_exp=(GCSs_in_TADs+GCSs_in_interTADs)*prob_of_interTAD_border
            #Bad TAD defenition
            GCSs_in_bad_TADs=TAD_assoc_dict[a]['GCSs in bad TAD']
            
            #Writing the output
            fileout.write('GCSs association with TAD and interTAD regions\n' + a + '\t' + TAD_inpath + '\n')
            print('Number of ' + str(a) + ' GCSs: ' + str(len(gcss_set)))
            print('Number of GCSs in TADs: ' + str(GCSs_in_TADs))
            print('Number of GCSs in interTADs: ' + str(GCSs_in_interTADs))
            print('Number of GCSs in TAD borders: ' + str(GCSs_in_TADs_border))
            print('Number of GCSs in interTAD borders: ' + str(GCSs_in_interTAD_border))
            print('Number of GCSs in bad TADs: ' + str(GCSs_in_bad_TADs))
            print('Corrected len of the genome (on deletions and on bad TADs): ' + str(TAD_genome_len))
            print('Len of TADs: ' + str(sum_TAD_len))
            print('Len of interTADs: ' + str(sum_interTAD_len))
            print('Binom test p-value for the number of GCSs in TADs: ')
            print(binom.cdf(GCSs_in_TADs, GCSs_in_TADs+GCSs_in_interTADs, prob_of_TAD))            
            print('Binom test p-value for the number of GCSs in interTADs: ')
            print(binom.cdf(GCSs_in_interTADs, GCSs_in_TADs+GCSs_in_interTADs, prob_of_interTAD))
            print('Binom test p-value for the number of GCSs in TAD borders: ')
            print(binom.cdf(GCSs_in_TADs_border, GCSs_in_TADs+GCSs_in_interTADs, prob_of_TAD_border))
            print('Binom test p-value for the number of GCSs in interTAD borders: ')
            print(binom.cdf(GCSs_in_interTAD_border, GCSs_in_TADs+GCSs_in_interTADs, prob_of_interTAD_border))            
            print('\n')
            fileout.write('Number of ' + str(a) + ' GCSs\t' + str(len(gcss_set)) + 
                          '\nWidth of the border region (nt)\t' + str(win_width) + 
                          '\nNumber of GCSs in TADs\t' + str(GCSs_in_TADs) + 
                          '\nNumber of GCSs in interTADs\t' + str(GCSs_in_interTADs) + 
                          '\nNumber of GCSs in TAD borders\t' + str(GCSs_in_TADs_border) +
                          '\nNumber of GCSs in interTAD borders\t' + str(GCSs_in_interTAD_border) + 
                          '\nNumber of GCSs in bad TADs\t' + str(GCSs_in_bad_TADs) + 
                          '\nCorrected len of the genome (on deletions and on bad TADs)\t' + str(TAD_genome_len) + 
                          '\nLen of TADs\t' + str(sum_TAD_len) + '\nLen of interTADs\t' + str(sum_interTAD_len) + 
                          '\nBinom test p-value for the number of GCSs in TADs\t' + str(binom.cdf(GCSs_in_TADs, GCSs_in_TADs+GCSs_in_interTADs, prob_of_TAD)) + 
                          '\nBinom test p-value for the number of GCSs in interTADs\t' + str(binom.cdf(GCSs_in_interTADs, GCSs_in_TADs+GCSs_in_interTADs, prob_of_interTAD)) + 
                          '\nChi-square test for the GCSs number asimmetry between TADs and interTADs\t' + str(chisquare([GCSs_in_TADs, GCSs_in_interTADs], f_exp=[GCSs_in_TADs_exp, GCSs_in_interTADs_exp])) +
                          '\nBinom test p-value for the number of GCSs in TAD borders\t' + str(binom.cdf(GCSs_in_TADs_border, GCSs_in_TADs+GCSs_in_interTADs, prob_of_TAD_border)) +
                          '\nBinom test p-value for the number of GCSs in interTAD borders\t' + str(binom.cdf(GCSs_in_interTAD_border, GCSs_in_TADs+GCSs_in_interTADs, prob_of_interTAD_border)) +
                          '\nBinom test p-value for the number of GCSs in TAD/interTAD borders\t' + str(binom.cdf(GCSs_in_TADs_border+GCSs_in_interTAD_border, GCSs_in_TADs+GCSs_in_interTADs, prob_of_TAD_border+prob_of_interTAD_border)) +
                          '\nChi-square test for the GCSs number asimmetry between TADs and interTADs borders\t' + str(chisquare([GCSs_in_TADs_border, GCSs_in_interTAD_border], f_exp=[GCSs_in_TADs_border_exp, GCSs_in_interTADs_border_exp])) +
                          '\n\n')
        TAD_score=[]
        InterTAD_score=[]
        for i in range(len(data_to_coloc)): #Iterates TADs
            TAD_score+=transcription[data_to_coloc[i][0]:data_to_coloc[i][1]]
            status=0
            for bad_TAD in Bad_TAD_reg: #Iterates bad TADs
                if bad_TAD[1]>=data_to_coloc[i][1]+10>=bad_TAD[0]: #bad TAD interval
                    status=1
                    break
            if status==0 and i<len(data_to_coloc)-1:
                InterTAD_score+=transcription[data_to_coloc[i][1]:data_to_coloc[i+1][0]]
            elif status==0 and i==len(data_to_coloc)-1:
                InterTAD_score+=transcription[data_to_coloc[i][1]:4647454]
        fileout.write('Comparison of TAD and InterTAD mean transcription level\n')
        print('Len of TAD regions: ' + str(len(TAD_score)))
        print('Len of interTAD regions: ' + str(len(InterTAD_score)))
        print('Mean TAD transcription level: ' + str(np.mean(TAD_score)))
        print('Mean interTAD transcription level: ' + str(np.mean(InterTAD_score)))
        transcription_t_test=stats.ttest_ind(TAD_score, InterTAD_score)
        print(transcription_t_test)
        print('\n')
        fileout.write('Mean TAD transcription level\t' + str(np.mean(TAD_score)) + 
                      '\nMean interTAD transcription level\t' + str(np.mean(InterTAD_score)) +
                      '\nt-test statistic\t' + str(transcription_t_test[0]) + 
                      '\nt-test p-value\t' + str(transcription_t_test[1]) + 
                      '\n\n')
    
    #Correction for deletions
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]]
    for i in deletions:
        genome_len+=i[0]-i[1] #Genome len correction on deletions

    if data_type=='Mut':
        Mut_assoc_dict={}
        for a, gcss_set in GCSs_sets_dict.items(): #Itarates GCSs sets (conditions)
            Mut_assoc_dict[a]={'Precise match':{}, '4bp GCSs':{}, '5bp vicinity':{}}
            for k, v in gcss_set.items(): #Itarates GCSs
                for Mut, mut_info in data_to_coloc.items(): #Iterates mutations
                    if k==Mut: #Mutation matches precisely GCS
                        Mut_assoc_dict[a]['Precise match'][k]=[v, mut_info]
                    elif k+4>Mut>=k: #Mutation falls into 4bp GCS region
                        Mut_assoc_dict[a]['4bp GCSs'][k]=[v, mut_info]
                    elif k+9>=Mut>=k-5: #Mutation falls into 5bp vicinity of GCSs 4bp region
                        Mut_assoc_dict[a]['5bp vicinity'][k]=[v, mut_info]
            fileout.write('GCSs association with mutations\n' + a + '\t' + Mut_inpath + '\n')
            print('Number of ' + str(a) + ' GCSs: ' + str(len(gcss_set)))
            print('Corrected len of the genome (on deletions): ' + str(genome_len))
            print('Total number of mutations: ' + str(len(data_to_coloc)))
            print('Number of mutations precisely match GCSs: ' + str(len(Mut_assoc_dict[a]['Precise match'])))
            print('Binom test p-value for the number of mutations that precisely match GCSs: ')
            print(binom.cdf(len(Mut_assoc_dict[a]['Precise match']), len(data_to_coloc), len(gcss_set)/genome_len))
            print('Number of mutations that fall into 4bp GCSs regions: ' + str(len(Mut_assoc_dict[a]['4bp GCSs'])))
            print('Binom test p-value for the number of mutations that fall into 4bp GCSs regions: ')
            print(binom.cdf(len(Mut_assoc_dict[a]['4bp GCSs']), len(data_to_coloc), len(gcss_set)*4/genome_len))
            print('Number of mutations that fall into 5bp vicinity of 4bp GCSs region: ' + str(len(Mut_assoc_dict[a]['5bp vicinity'])))
            print('Binom test p-value for the number of mutations that fall into 5bp vicinity of 4bp GCSs region: ')
            print(binom.cdf(len(Mut_assoc_dict[a]['5bp vicinity']), len(data_to_coloc), len(gcss_set)*14/genome_len))            
            print('\n')
            fileout.write('Number of ' + str(a) + ' GCSs\t' + str(len(gcss_set)) + '\nCorrected len of the genome (on deletions)\t' + str(genome_len) + 
                          '\nTotal number of mutations\t' + str(len(data_to_coloc)) + 
                          '\nNumber of mutations precisely match GCSs\t' + str(len(Mut_assoc_dict[a]['Precise match'])) + 
                          '\nBinom test p-value for the number of mutations that precisely match GCSs\t' + str(binom.cdf(len(Mut_assoc_dict[a]['Precise match']), len(data_to_coloc), len(gcss_set)/genome_len)) + 
                          '\nNumber of mutations that fall into 4bp GCSs regions\t' + str(len(Mut_assoc_dict[a]['4bp GCSs'])) + 
                          '\nBinom test p-value for the number of mutations that fall into 4bp GCSs regions\t' + str((binom.cdf(len(Mut_assoc_dict[a]['4bp GCSs']), len(data_to_coloc), len(gcss_set)*4/genome_len))) + 
                          '\nNumber of mutations that fall into 5bp vicinity of 4bp GCSs regions\t' + str(len(Mut_assoc_dict[a]['5bp vicinity'])) +
                          '\nBinom test p-value for the number of mutations that fall into 5bp vicinity of 4bp GCSs regions\t' + str(binom.cdf(len(Mut_assoc_dict[a]['5bp vicinity']), len(data_to_coloc), len(gcss_set)*14/genome_len)) + 
                          '\n\n')
    fileout.close()
    return

#######
#Funcrions wrapper.
#######

def wrapper(input_dict, transcription_path, TAD_inpath, Mut_inpath, Bad_TAD_reg, outpath):
    #Data parsing
    GCSs_sets_dict=trusted_GCSs_parsing(input_dict)
    transcription=transcription_data_parser(transcription_path)
    TAD_coord=TAD_parser(TAD_inpath)
    Mut_dict=Mutations_parser(Mut_inpath)
    #Data classification
    colocal_search(GCSs_sets_dict, transcription, 'TAD', TAD_coord, Bad_TAD_reg, TAD_inpath, Mut_inpath, outpath)
    colocal_search(GCSs_sets_dict, transcription, 'Mut', Mut_dict, Bad_TAD_reg, TAD_inpath, Mut_inpath, outpath)
    return

Bad_TAD_reg_150kb=[[220002, 394224], [776546, 966544], [1147701, 1307836], [3420979, 3432802], [4196538, 4315800]] #Contains deletions
Bad_TAD_reg_80kb=[[220002, 379224], [776546, 876544], [1177701, 1307836], [3420979, 3432802], [4196538, 4279712]] #Contains deletions
wrapper(path_to_GCSs_files, Transcription_data_path, TAD_150_inputpath, Mut_wt_inputpath, Bad_TAD_reg_150kb, Outputpath)
wrapper(path_to_GCSs_files, Transcription_data_path, TAD_80_inputpath, Mut_mut_inputpath, Bad_TAD_reg_80kb, Outputpath)

print('Script ended its work succesfully!') 