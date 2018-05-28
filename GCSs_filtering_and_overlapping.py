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
path_to_cfx_replicas={'Cfx_1': '',
                      'Cfx_2': '',
                      'Cfx_3': ''}
path_to_rifcfx_replicas={'RifCfx_1': '',
                      'RifCfx_2': '',
                      'RifCfx_3': ''}
path_to_microcin_replicas={'Micro_1': '',
                      'Micro_2': '',
                      'Micro_3': ''}
path_to_oxo_replicas={'Oxo_1': '',
                      'Oxo_2': '',
                      'Oxo_3': ''}
#Configuration of the output for the GCSs data in replicas.
Replicas_path_out=''
Cfx_name='Cfx'
RifCfx_name='RifCfx'
Micro_name='Micro'
Oxo_name='Oxo'
#Configuration of the output for GCSs trusted.
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
#Parsing raw GCSs coordinates, returns dictionary - GCSs_coordinate:N3E.
#######

def read_GCSs_file(GCSs_file_path):
    GCSs_dict={}
    GCSs_in=open(GCSs_file_path, 'r')
    for line in GCSs_in:
        line=line.rstrip().split('\t')
        if line[0] not in ['GCSs_coordinate']:
            GCSs_dict[int(line[0])]=float(line[1])
    GCSs_in.close()
    return GCSs_dict

#######
#Combines replicates into one GCSs table.
#######

def combine_replicates(replicas_dict, path_out, name):
    #Merges a range of replicates
    GCSs_replicas_dict={}
    names_ar=[]
    for key, value in replicas_dict.items():
        names_ar.append(k)
        for k, v in read_GCSs_file(value).items():
            if k in GCSs_replicas_dict:
                GCSs_replicas_dict[k].append(v)
            else:
                add_el=[]
                for j in range(i):
                    add_el.append(0)
                add_el.append(v)
                GCSs_replicas_dict[k]=add_el
        for k, v in GCSs_replicas_dict.items():
            if len(v)<i+1:
                GCSs_replicas_dict[k].append(0)
    #Writes merged GCSs data
    fileout=open(path_out + name + '_GCSs_replicates.txt', 'w')
    fileout.write('GCSs_coordinate\t')
    for i in names_ar:
        fileout.write(str(i) + '_N3E\t')
    fileout.write('\n')
    for k, v in GCSs_replicas_dict.items():
        fileout.write(str(k) + '\t')
        for i in v:
            fileout.write(str(v) + '\t')
        fileout.write('\n')
    fileout.close()
    return GCSs_replicas_dict
        

#######
#Returns only trusted GCSs - observed at least 2 times within 3 biological replicates.
#Data organization: 1. coordinate of GCSs, 2.-4. N3E values for biological replicates 1-3
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

def trusted_GCSs_calling(GCSs_dictionary):
    ar=[]
    for k, v in GCSs_dictionary.items():
        if trusted(v)!="No signal":
            ar.append([k, trusted(v)])
    return ar

def replicas_comb_trust_wrapper(replicas_dict, path_out, name):
    print('Now working with: ' + str(name))
    cur_GCSs_dict=combine_replicates(replicas_dict, path_out, name)
    cur_GCSs_trusted=trusted_GCSs_calling(cur_GCSs_dict)
    print('Number of trusted GCSs for ' + str(name) + ' : ' + str(len(cur_GCSs_trusted)))
    return cur_GCSs_trusted

Cfx=replicas_comb_trust_wrapper(path_to_cfx_replicas, Replicas_path_out, Cfx_name)
RifCfx=replicas_comb_trust_wrapper(path_to_rifcfx_replicas, Replicas_path_out, RifCfx_name)
Micro=replicas_comb_trust_wrapper(path_to_microcin_replicas, Replicas_path_out, Micro_name)
Oxo=replicas_comb_trust_wrapper(path_to_oxo_replicas, Replicas_path_out, Oxo_name)

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
#GCSs sets average N3E estimation.
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
#Write down files with GCSs lists - trusted or shared.
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
 
print('Script ended its work succesfully!') 
 
