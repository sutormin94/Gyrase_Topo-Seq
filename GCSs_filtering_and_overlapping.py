###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes raw GCSs data, returns only trusted GCSs, 
#computes GCSs shared between different conditions, 
#draws Venn diagrams of the sets overlappings, 
#writes GCSs sets.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import collections
from matplotlib_venn import venn2, venn3, venn3_circles
import numpy as np

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Path to the working directory
pwd="/data/Gyrase/Data_preparation/"

#Input data
path_to_cfx_replicas={'Cfx_1': pwd + "Cfx_10mkM/GCSs_calling_0_05/Cfx_1/Cfx_1_raw_GCSs_called.txt",
                      'Cfx_2': pwd + "Cfx_10mkM/GCSs_calling_0_05/Cfx_2/Cfx_2_raw_GCSs_called.txt",
                      'Cfx_3': pwd + "Cfx_10mkM/GCSs_calling_0_05/Cfx_3/Cfx_3_raw_GCSs_called.txt"}
path_to_rifcfx_replicas={'RifCfx_1': pwd + "RifCfx/GCSs_calling_0_05_Un/RifCfx_1/RifCfx_1_raw_GCSs_called.txt",
                         'RifCfx_2':  pwd + "RifCfx/GCSs_calling_0_05_Un/RifCfx_2/RifCfx_2_raw_GCSs_called.txt",
                      'RifCfx_3':  pwd + "RifCfx/GCSs_calling_0_05_Un/RifCfx_3/RifCfx_3_raw_GCSs_called.txt"}
path_to_microcin_replicas={'Micro_1':  pwd + "Micro/GCSs_calling_0_05/Micro_1/Micro_1_raw_GCSs_called.txt",
                           'Micro_2': pwd + "Micro/GCSs_calling_0_05/Micro_2/Micro_2_raw_GCSs_called.txt",
                      'Micro_3': pwd + "Micro/GCSs_calling_0_05/Micro_3/Micro_3_raw_GCSs_called.txt"}
path_to_oxo_replicas={'Oxo_1': pwd + "Oxo/GCSs_calling_0_05/Oxo_1/Oxo_1_raw_GCSs_called.txt",
                      'Oxo_2': pwd + "Oxo/GCSs_calling_0_05/Oxo_2/Oxo_2_raw_GCSs_called.txt",
                      'Oxo_3': pwd + "Oxo/GCSs_calling_0_05/Oxo_3/Oxo_3_raw_GCSs_called.txt"}
#Input data in one dict (for one output table contains all replices and raw GCSs)
path_to_replicas={'Cfx_1': pwd + "Cfx_10mkM/GCSs_calling_0_05/Cfx_1/Cfx_1_raw_GCSs_called.txt",
                  'Cfx_2': pwd + "Cfx_10mkM/GCSs_calling_0_05/Cfx_2/Cfx_2_raw_GCSs_called.txt",
                      'Cfx_3': pwd + "Cfx_10mkM/GCSs_calling_0_05/Cfx_3/Cfx_3_raw_GCSs_called.txt",
                      'RifCfx_1': pwd + "RifCfx/GCSs_calling_0_05_Un/RifCfx_1/RifCfx_1_raw_GCSs_called.txt",
                      'RifCfx_2':  pwd + "RifCfx/GCSs_calling_0_05_Un/RifCfx_2/RifCfx_2_raw_GCSs_called.txt",
                      'RifCfx_3':  pwd + "RifCfx/GCSs_calling_0_05_Un/RifCfx_3/RifCfx_3_raw_GCSs_called.txt",
                      'Micro_1':  pwd + "Micro/GCSs_calling_0_05/Micro_1/Micro_1_raw_GCSs_called.txt",
                      'Micro_2': pwd + "Micro/GCSs_calling_0_05/Micro_2/Micro_2_raw_GCSs_called.txt",
                      'Micro_3': pwd + "Micro/GCSs_calling_0_05/Micro_3/Micro_3_raw_GCSs_called.txt",
                      'Oxo_1': pwd + "Oxo/GCSs_calling_0_05/Oxo_1/Oxo_1_raw_GCSs_called.txt",
                      'Oxo_2': pwd + "Oxo/GCSs_calling_0_05/Oxo_2/Oxo_2_raw_GCSs_called.txt",
                      'Oxo_3': pwd + "Oxo/GCSs_calling_0_05/Oxo_3/Oxo_3_raw_GCSs_called.txt"}

#Configuration of the output for the GCSs data in replicas.
Replicas_path_out="/data/Gyrase/GCSs_sets/"
if not os.path.exists(Replicas_path_out):
    os.makedirs(Replicas_path_out)
Cfx_name="Cfx_10mkM"
RifCfx_name="RifCfx"
Micro_name="Micro"
Oxo_name="Oxo"
All_conditions_name="All_conditions_GCSs"
#Configuration of the output for GCSs trusted.
Cfx_path=Replicas_path_out + "Cfx_10mkM_trusted_GCSs.txt"
RifCfx_path=Replicas_path_out + "RifCfx_trusted_GCSs.txt"
Micro_path=Replicas_path_out + "Micro_trusted_GCSs.txt"
Oxo_path=Replicas_path_out + "Oxo_trusted_GCSs.txt"
Cfx_Micro_path=Replicas_path_out + "Cfx_10mkM_Micro_shared_trusted_GCSs.txt"
Cfx_Oxo_path=Replicas_path_out + "Cfx_10mkM_Oxo_shared_trusted_GCSs.txt"
Micro_Oxo_path=Replicas_path_out + "Micro_Oxo_shared_trusted_GCSs.txt"
Cfx_Micro_Oxo_path=Replicas_path_out + "Cfx_10mkM_Micro_Oxo_shared_trusted_GCSs.txt"
Cfx_RifCfx_shared_GCSs_path=Replicas_path_out + "Cfx_10mkM_RifCfx_shared_trusted_GCSs.txt"
#Outpath for Venn diagrams.
plot_outpath=Replicas_path_out

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
    for key, value in replicas_dict.items(): #Iterates replicas
        names_ar.append(key)
        #Read file with raw GCSs
        Raw_GCSs_dict=read_GCSs_file(value)
        for k, v in Raw_GCSs_dict.items(): #Iterates raw GCSs
            #Table filling process initiation
            if len(names_ar)==1:
                GCSs_replicas_dict[k]=[v]
            #Table filling process continuing (the table already contains at least one GCSs set)
            else:
                #If GCSs is already in the table
                if k in GCSs_replicas_dict:
                    GCSs_replicas_dict[k].append(v)
                #If this is the first occurrence of the element in a NON empty table.
                else:
                    add_el=[]
                    for j in range(len(names_ar)-1):
                        add_el.append(0)
                    add_el.append(v)
                    GCSs_replicas_dict[k]=add_el
        #If table body line contains less elements than header does, hence add zero.
        for k, v in GCSs_replicas_dict.items():
            if len(v)<len(names_ar):
                GCSs_replicas_dict[k].append(0)
    #Sorting the list of dictionary keys.
    GCSs_replicas_dict_sorted=collections.OrderedDict(sorted(GCSs_replicas_dict.items()))
    #Writes merged GCSs data
    fileout=open(path_out + name + '_GCSs_replicates.txt', 'w')
    #Header
    fileout.write('GCSs_coordinate\t')
    for i in names_ar:
        fileout.write(str(i) + '_N3E\t')
    fileout.write('\n')
    #Body of the table
    for k, v in GCSs_replicas_dict_sorted.items():
        fileout.write(str(k) + '\t')
        for i in GCSs_replicas_dict_sorted[k]:
            fileout.write(str(i) + '\t')
        fileout.write('\n')
    fileout.close()
    return GCSs_replicas_dict
 
#Prepares GCSs table for all conditions
combine_replicates(path_to_replicas, Replicas_path_out, All_conditions_name)      

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
#GCSs shared between pairs of antibiotics - Cfx, Micro and Oxo and between Cfx and RifCfx.
#######

def pairs_construction(ar1, ar2):
    double=[]
    for i in range(len(ar1)):
        for j in range(len(ar2)):
            if ar1[i][0]==ar2[j][0]:
                double.append([ar1[i][0], ar1[i][1], ar2[j][1]]) #GCSs coordinate, N3E_1, N3E_2 
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
                triple.append([ar12[i][0], ar12[i][1], ar12[i][2], ar3[j][1]]) #GCSs coordinate, N3E_1, N3E_2, N3E_3
    return triple

Cfx_Micro_Oxo_shared_GCSs=triple_construction(Cfx_Micro_shared_GCSs, Oxo)
print('Number of GCSs shared between Cfx, Micro and Oxo: ' + str(len(Cfx_Micro_Oxo_shared_GCSs)) +'\n')

#######
#Parses replicas, overlaps lists of GCSs, output data for Venn diagram construction.
#######

def replicates_parsing_to_list_and_overlapping(replicas_dict, name):
    #Parsing
    GCSs_dict={}
    for k, v in replicas_dict.items(): #Iterate replicas.
        GCSs_dict[k]=[]
        for c, h in read_GCSs_file(v).items(): #Iterate GCSs.
            GCSs_dict[k].append([c, h])
    #Overlapping
    one_two=pairs_construction(GCSs_dict[name+str(1)], GCSs_dict[name+str(2)])
    one_three=pairs_construction(GCSs_dict[name+str(1)], GCSs_dict[name+str(3)])
    two_three=pairs_construction(GCSs_dict[name+str(2)], GCSs_dict[name+str(3)])
    one_two_three=triple_construction(one_two, GCSs_dict[name+str(3)])
    #Venn input description (for 3 sets): one, two, three, one_two, one_three, two_three, one_two_three
    venn_input=[len(GCSs_dict[name+str(1)])-len(one_two)-len(one_three)+len(one_two_three), 
                len(GCSs_dict[name+str(2)])-len(one_two)-len(two_three)+len(one_two_three), 
                len(GCSs_dict[name+str(3)])-len(one_three)-len(two_three)+len(one_two_three),
                len(one_two)-len(one_two_three), len(one_three)-len(one_two_three), len(two_three)-len(one_two_three),
                len(one_two_three)]
    return venn_input

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

venn2(subsets = (venn_data_2), set_labels = ("Ciprofloxacin", "Rifampicin Ciprofloxacin"))
plt.savefig(plot_outpath+'Cfx_RifCfx_venn.png', dpi=320)
plt.close()

print("Cfx Micro Oxo subsets volumes: " + str(venn_data_3))
venn3(subsets = (venn_data_3), set_labels = ('Ciprofloxacin', 'Microcin B17', 'Oxolinic acid'))
plt.savefig(plot_outpath+'Cfx_Micro_Oxo_venn.png', dpi=320)
plt.close()

venn3(subsets = (replicates_parsing_to_list_and_overlapping(path_to_cfx_replicas, 'Cfx_')), set_labels = ('Cfx_1', 'Cfx_2', 'Cfx_3'))
plt.savefig(plot_outpath+'Cfx_replicas_venn.png', dpi=320)
plt.close()

venn3(subsets = (replicates_parsing_to_list_and_overlapping(path_to_rifcfx_replicas, 'RifCfx_')), set_labels = ('RifCfx_1', 'RifCfx_2', 'RifCfx_3'))
plt.savefig(plot_outpath+'RifCfx_replicas_venn.png', dpi=320)
plt.close()

venn3(subsets = (replicates_parsing_to_list_and_overlapping(path_to_microcin_replicas, 'Micro_')), set_labels = ('Micro_1', 'Micro_2', 'Micro_3'))
plt.savefig(plot_outpath+'Micro_replicas_venn.png', dpi=320)
plt.close()

venn3(subsets = (replicates_parsing_to_list_and_overlapping(path_to_oxo_replicas, 'Oxo_')), set_labels = ('Oxo_1', 'Oxo_2', 'Oxo_3'))
plt.savefig(plot_outpath+'Oxo_replicas_venn.png', dpi=320)
plt.close()

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
               Micro_Oxo_path: Antibs_GCSs_sets_pair_shared[2],
               Cfx_Micro_Oxo_path: Cfx_Micro_Oxo_shared_GCSs}

def write_GCSs_file(dictionary):
    for k, v in dictionary.items(): #Iterates lists to be written
        v.sort(key=lambda tup: tup[0])  #Sorting lists by the zero elements of the sublists they consist of 
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
    ar.sort(key=lambda tup: tup[0])
    for i in range(len(ar)):
        fileout.write(str(ar[i][0]) + '\t' + str(ar[i][1]) + '\t' + str(ar[i][2]) + '\n')
    fileout.close()
    return
    
write_Cfx_RifCfx_shared_GCSs(Cfx_RifCfx_shared_GCSs, Cfx_RifCfx_shared_GCSs_path)
 
print('Script ended its work succesfully!') 
 
