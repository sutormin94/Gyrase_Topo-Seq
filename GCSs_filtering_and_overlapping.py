###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes raw GCSs data, returns only trusted GCSs, 
#computes GCSs shared between different conditions, 
#draws Venn diagrams of the sets overlappings, 
#writes GCSs sets.
###############################################

#######
#Packages to be imported
#######

import matplotlib.pyplot as plt
import numpy as np

#######
#Variables to be defined
#######

print('Variables to be defined:')

#Input data
path_to_cfx_raw_data=''
path_to_rifcfx_raw_data=''
path_to_microcin_raw_data=''
path_to_oxo_raw_data=''
#Output_data
Cfx_path=''
RifCfx_path=''
Micro_path=''
Oxo_path=''
Cfx_Micro_path=''
Cfx_Oxo_path=''
Micro_Oxo_path=''
Cfx_Micro_Oxo_path=''
Cfx_RifCfx_shared_GCSs_path=''


#######
#Parsing raw GCSs coordinates, returns only trusted GCSs - observed at least 2 times within 3 biological replicates.
#Raw data is organized as a tabbed table, consist of 4 columns: 1. coordinate of GCSs, 2.-4. N3E values for biological replicates 1-3
#######

def trusted(ar):
    av_height=0
    ind=0
    for i in range(len(ar)):
        if ar[i]>0:
            ind=ind+1
            av_height=av_height+ar[i]
    if ind>1:
        return av_height/ind
    else:
        return "No signal"

def trusted_GCSs_calling(path):
    ar=[]
    filein=open(path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ["peak"]:
            line=[int(line[0]), float(line[1]), float(line[2]), float(line[3])]
        if trusted(line[1:])!="No signal":
            ar.append([line[0], trusted(line[1:])])
    filein.close()
    return ar

Cfx=trusted_GCSs_calling(path_to_cfx_raw_data)
RifCfx=trusted_GCSs_calling(path_to_rifcfx_raw_data)
Micro=trusted_GCSs_calling(path_to_microcin_raw_data)
Oxo=trusted_GCSs_calling(path_to_oxo_raw_data)

print('Number of trusted Cfx GCSs: ' + str(len(Cfx)))
print('Number of trusted RifCfx GCSs: ' + str(len(RifCfx)))
print('Number of trusted Micro GCSs: ' + str(len(Micro)))
print('Number of trusted Oxo GCSs: ' + str(len(Oxo)) + '\n')

Antibs_GCSs_sets=[Cfx, RifCfx, Micro, Oxo]

#######
#GCSs shared between pairs of antibiotics - Cfx, Micro and Oxo and between Cfx and RifCfx
#######

def pairs_construction(ar1, ar2):
    double=[]
    for i in range(len(ar1)):
        for j in range(len(ar2)):
            if ar1[i][0]==ar2[j][0]:
                cp_info=[]
                cp_info.append(ar1[i][0])
                cp_info.append(ar1[i][1])
                cp_info.append(ar2[j][1])
                double.append(cp_info)      
    return double

Cfx_RifCfx_shared_GCSs=pairs_construction(Cfx, RifCfx)
print('Number of GCSs shared between Cfx and RifCfx: ' + str(len(Cfx_RifCfx_shared_GCSs)) + '\n')

Cfx_Micro_shared_GCSs=pairs_construction(Cfx, Micro)
Cfx_Oxo_shared_GCSs=pairs_construction(Cfx, Oxo)
Micro_Oxo_shared_GCSs=pairs_construction(Micro, Oxo)

print('Number of GCSs shared between Cfx and Micro: ' + str(len(Cfx_Micro_shared_GCSs)))
print('Number of GCSs shared between Cfx and Oxo: ' + str(len(Cfx_Oxo_shared_GCSs)))
print('Number of GCSs shared between Micro and Oxo: ' + str(len(Micro_Oxo_shared_GCSs)) + '\n')

Antibs_GCSs_sets_pair_shared=[Cfx_Micro_shared_GCSs, Cfx_Oxo_shared_GCSs, Micro_Oxo_shared_GCSs]

#######
#GCSs shared between 3 antibiotics
#######

def triple_construction(ar12, ar3):
    triple=[]
    for i in range(len(ar12)):
        for j in range(len(ar3)):
            if ar12[i][0]==ar3[j][0]:
                ct_info=[]
                ct_info.append(ar12[i][0])
                ct_info.append(ar12[i][1])
                ct_info.append(ar12[i][2])
                ct_info.append(ar3[j][1])
                triple.append(ct_info)      
    return triple

Cfx_Micro_Oxo_shared_GCSs=triple_construction(Cfx_Micro_shared_GCSs, Oxo)
print('Number of GCSs shared between Cfx, Micro and Oxo: ' + str(len(Cfx_Micro_Oxo_shared_GCSs)) +'\n')

#######
#Venn diagram represents GCSs sets overlapping.
#description2: one, two, one_two
#description3: one, two, three, one_two, one_three, two_three, one_two_three
#######

venn_data_2=[len(Cfx)-len(Cfx_RifCfx_shared_GCSs), len(RifCfx)-len(Cfx_RifCfx_shared_GCSs), len(Cfx_RifCfx_shared_GCSs)]
venn_data_3=[len(Cfx)-len(Cfx_Micro_shared_GCSs)-len(Cfx_Oxo_shared_GCSs)+len(Cfx_Micro_Oxo_shared_GCSs), 
             len(Micro)-len(Cfx_Micro_shared_GCSs)-len(Micro_Oxo_shared_GCSs)+len(Cfx_Micro_Oxo_shared_GCSs), 
             len(Oxo)-len(Cfx_Oxo_shared_GCSs)-len(Micro_Oxo_shared_GCSs)+len(Cfx_Micro_Oxo_shared_GCSs),
             len(Cfx_Micro_shared_GCSs)-len(Cfx_Micro_Oxo_shared_GCSs), 
             len(Cfx_Oxo_shared_GCSs)-len(Cfx_Micro_Oxo_shared_GCSs), 
             len(Micro_Oxo_shared_GCSs)-len(Cfx_Micro_Oxo_shared_GCSs), 
             len(Cfx_Micro_Oxo_shared_GCSs)]

venn2(subsets = (venn_data_2[0], venn_data_2[1], venn_data_2[2]), set_labels = ("Ciprofloxacin", "Rifampicin Ciprofloxacin"))
plt.show()
venn3(subsets = (venn_data_3[0], venn_data_3[1], venn_data_3[2], venn_data_3[3], venn_data_3[4], venn_data_3[5], venn_data_3[6]), set_labels = ('Ciprofloxacin', 'Microcin B17', 'Oxolinic acid'))
plt.show()

#######
#GCSs sets average N3E estimation
#######

def average_height(ar):
    av_he=0
    for i in range(len(ar)):
        peak_he=np.mean(ar[i][1:])
        av_he=av_he+peak_he
    return av_he/len(ar)

print('Cfx average GCSs N3E: ' + str(average_height(Cfx)))
print('Micro average GCSs N3E: ' + str(average_height(Micro)))
print('Oxo average GCSs N3E: ' + str(average_height(Oxo)))
print('Cfx and Micro average GCSs N3E: ' + str(average_height(Cfx_Micro_shared_GCSs)))
print('Cfx and Oxo average GCSs N3E: ' + str(average_height(Cfx_Oxo_shared_GCSs)))
print('Micro and Oxo average GCSs N3E: ' + str(average_height(Micro_Oxo_shared_GCSs)))
print('Cfx, Micro and Oxo average GCSs N3E: ' + str(average_height(Cfx_Micro_Oxo_shared_GCSs)) + '\n')


#######
#Write down files with GCSs lists - trusted or shared
#######

All_GCSs_sets={Cfx_path: Antibs_GCSs_sets[0],
               RifCfx_path: Antibs_GCSs_sets[1],
               Micro_path: Antibs_GCSs_sets[2],
               Oxo_path: Antibs_GCSs_sets[3],
               Cfx_Micro_path: Antibs_GCSs_sets_pair_shared[0],
               Cfx_Oxo_path: Antibs_GCSs_sets_pair_shared[1],
               Micro_Oxo_path: Antibs_GCSs_sets_pair_shared[3],
               Cfx_Micro_Oxo_path: Cfx_Micro_Oxo_shared_GCSs}

def write_GCSs_file(dictionary):
    for k, v in dictionary.items():
        fileout=open(k, 'w')
        fileout.write('GCSs_coordinate\tN3E\n')
        for i in range(len(v)):
            fileout.write(str(v[i][0]) + '\t' + str(np.mean(v[i][1:])) + '\n')
        fileout.close()
    return

write_GCSs_file(All_GCSs_sets)


def write_Cfx_RifCfx_shared_GCSs(ar, path):
    fileout=open(path, 'w')
    fileout.write('GCSs_coordinate\tCfx_N3E\tRifCfx_N3E\n')
    for i in range(len(ar)):
        fileout.write(str(ar[i][0]) + '\t' + str(ar[i][1]) + '\t' + str(ar[i][2]) + '\n')
    fileout.close()
    return
    
write_Cfx_RifCfx_shared_GCSs(Cfx_RifCfx_shared_GCSs, Cfx_RifCfx_shared_GCSs_path)
 
print('Script ended succesfully!') 
 
