###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script analysis colocalization of mutations (Foster, 2015) and TADs borders (Lioy, 2018) with GCSs.
###############################################

#######
#Packages to be imported.
#######

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import binom

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Input data - GCSs, TAB.
path_to_GCSs_files={'Cfx': 'C:\Sutor\science\DNA-gyrase\Results\Final_data_2\GCSs_sets\Cfx_trusted_one_position.txt',
                    'RifCfx': 'C:\Sutor\science\DNA-gyrase\Results\Final_data_2\GCSs_sets\Rif_Cfx_trusted_one_position.txt',
                    'Micro': 'C:\Sutor\science\DNA-gyrase\Results\Final_data_2\GCSs_sets\Micro_trusted_one_position.txt',
                    'Oxo': 'C:\Sutor\science\DNA-gyrase\Results\Final_data_2\GCSs_sets\Oxo_trusted_one_position.txt'}
#TAD coordinates.
TAD_150_inputpath='C:\Sutor\science\DNA-gyrase\Results\Hi-C\TAD_calling\CIDs_150kb_31_W3110_Mu_compatible.bed'
TAD_80_inputpath='C:\Sutor\science\DNA-gyrase\Results\Hi-C\TAD_calling\CIDs_80kb_57_W3110_Mu_compatible.bed'
#Mutations coordinates.
Mut_wt_inputpath='C:\Sutor\science\E_coli_variome\Foster_data\Foster,2015,WT_E_coli_w3110_Mu_compatible.vcf'
Mut_mut_inputpath='C:\Sutor\science\E_coli_variome\Foster_data\Foster,2015,Mut_E_coli_w3110_Mu_compatible.vcf'


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
#Colocalization.
#######

def colocal_search(GCSs_sets_dict, data_type, data_to_coloc, Bad_TAD_reg):
    genome_len=4647454 #Overall genome len
    deletions=[[274500, 372148], [793800, 807500], [1199000, 1214000]]
    for i in deletions:
        genome_len+=i[0]-i[1] #Genome len corrected on deletions
        
    if data_type=='TAD':
       
        #Genome len correction on bad TAD regions.
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
            TAD_assoc_dict[a]={'GCSs in TAD':0, 'GCSs in interTAD':0, 'GCSs in bad TAD':0}
            for k, v in gcss_set.items(): #Itarates GCSs
                status=0
                for TAD in data_to_coloc: #Iterates TADs
                    if TAD[1]>=k>=TAD[0]: #GCS falls into TAD
                        TAD_assoc_dict[a]['GCSs in TAD']+=1
                        status=1
                        break
                if status==0: #GCS not in TAD
                    for bad_TAD in Bad_TAD_reg:
                        if bad_TAD[1]>k>bad_TAD[0]: #GCSs fall in the bad TAD region.
                            status=1
                            TAD_assoc_dict[a]['GCSs in bad TAD']+=1
                            break
                    if status==0: #GCSs is in interTAD
                        TAD_assoc_dict[a]['GCSs in interTAD']+=1
            
            print('Number of ' + str(a) + ' GCSs: ' + str(len(gcss_set)))
            print('Number of GCSs in TADs: ' + str(TAD_assoc_dict[a]['GCSs in TAD']))
            print('Number of GCSs in interTADs: ' + str(TAD_assoc_dict[a]['GCSs in interTAD']))
            print('Number of GCSs in bad TADs: ' + str(TAD_assoc_dict[a]['GCSs in bad TAD']))
            print('Len of corrected genome: ' + str(TAD_genome_len))
            print('Len of TADs: ' + str(sum_TAD_len))
            print('Len of interTADs: ' + str(sum_interTAD_len) + '\n')
    
    elif data_type=='Mut':
        Mut_assoc_dict={}
        for a, gcss_set in GCSs_sets_dict.items(): #Itarates GCSs sets (conditions)
            Mut_assoc_dict[a]={}
            for k, v in gcss_set.items(): #Itarates GCSs
                for Mut, mut_info in data_to_coloc.items(): #Iterates mutations
                    if k==Mut:
                        Mut_assoc_dict[a][k]=[v, mut_info]
            print('Number of ' + str(a) + ' GCSs: ' + str(len(gcss_set)))
            print('Len of the genome: ' + str(genome_len))
            print('Total number of mutations: ' + str(len(data_to_coloc)))
            print('Number of mutations matched GCSs: ' + str(len(Mut_assoc_dict[a])) + '\n')
    return

#######
#Funcrions wrapper.
#######

def wrapper(input_dict, TAD_inpath, Mut_inpath, Bad_TAD_reg):
    #Data parsing
    GCSs_sets_dict=trusted_GCSs_parsing(input_dict)
    TAD_coord=TAD_parser(TAD_inpath)
    Mut_dict=Mutations_parser(Mut_inpath)
    #Data classification
    colocal_search(GCSs_sets_dict, 'TAD', TAD_coord, Bad_TAD_reg)
    colocal_search(GCSs_sets_dict, 'Mut', Mut_dict, Bad_TAD_reg)
    
    return

Bad_TAD_reg_150kb=[[1147701, 1307836], [3420979, 3432802], [4196538, 4315800]]
Bad_TAD_reg_80kb=[[1177701, 1307836], [3420979, 3432802], [4196538, 4279712]] 
wrapper(path_to_GCSs_files, TAD_150_inputpath, Mut_wt_inputpath, Bad_TAD_reg_150kb)
wrapper(path_to_GCSs_files, TAD_80_inputpath, Mut_mut_inputpath, Bad_TAD_reg_80kb)

print('Script ended its work succesfully!') 