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

#Input: REPs classification and relative positions within BIME-2 regions, TAB.
REPs_coord_input="C:\Sutor\science\DNA-gyrase\Results\Final_data_2\GCSs_association_with_REPs\REP_coordinates_within_BIMEs.txt"

#Input: GCSs, TAB.
path_to_GCSs_files={'Cfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\Cfx_10mkM_trusted_GCSs_h_s.txt",
                    'RifCfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\RifCfx_trusted_GCSs_h_s.txt",
                    'Micro': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\Micro_trusted_GCSs_h_s.txt",
                    'Oxo': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\Oxo_trusted_GCSs_h_s.txt"}  

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
#REPs coordinates parsing and REPs classification.
#######

def REPs_parsing(input_reps_path):
    reps_in=open(input_reps_path, "r")
    y_L=[]
    y_D=[]
    z2_L=[]
    z2_D=[]
    for line in reps_in:
        line=line.rstrip().split('\t')
        line[2]=line[2].split('#')
        for rep in line[2]:
            rep_info=rep.split(':')
            if rep_info[0]=='yL':
                y_L.append([int(line[0])+int(rep_info[1]), int(line[0])+int(rep_info[2])])
            if rep_info[0]=='yD':
                y_D.append([int(line[0])+int(rep_info[1]), int(line[0])+int(rep_info[2])])
            if rep_info[0]=='z2L':
                z2_L.append([int(line[0])+int(rep_info[1]), int(line[0])+int(rep_info[2])])
            if rep_info[0]=='z2D':
                z2_D.append([int(line[0])+int(rep_info[1]), int(line[0])+int(rep_info[2])])
    REPs_coords=[y_L, y_D, z2_L, z2_D]
    reps_in.close()
    return REPs_coords

#######
#GCSs localization within classified REPs.
#######

def locate_gcs(GCSs_sets_dict, REPs_coords):
    return

print('Script ended its work succesfully!') 