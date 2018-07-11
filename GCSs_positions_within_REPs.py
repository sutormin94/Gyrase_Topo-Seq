###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script analysis localization of GCSs within REP regions.
#Additionally it calculates mean scores for each position of REPs.
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
from Bio import SeqIO
from Bio.Seq import Seq

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

#Input: path to the E. coli genome (source of sequences for PFM/PWM construction), FASTA.
Genome_seq_path="C:\Sutor\science\DNA-gyrase\Genomes\E_coli_w3110_G_Mu.fasta"

#Input: path to the score track for the E. coli genome, WIG.
Score_data_path="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\Score_tracks\E_coli_w3110_G_Mu_score.wig"

#Output: localized GCSs, TAB.
Outputpath="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_within_REPs\\"
if not os.path.exists(Outputpath):
    os.makedirs(Outputpath)    
    
    
#######
#FASTA sequences parsing.
#######

def obtain_seq(seq_path):
    seq_oi=open(seq_path, 'r')
    for record in SeqIO.parse(seq_oi, "fasta"):
        sequence=str(record.seq)
    seq_oi.close()      
    return sequence, record.id

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
#REPs coordinates parsing and REPs classification.
#######

def REPs_parsing(input_reps_path, outpath):
    reps_in=open(input_reps_path, "r")
    rep_coord={'yL': [], 'yD':[], 'z2L':[], 'z2D':[]}
    rep_lens={'yL': [], 'yD':[], 'z2L':[], 'z2D':[]}
    for line in reps_in:
        line=line.rstrip().split('\t')
        line[2]=line[2].split('#')
        for rep in line[2]:
            rep_info=rep.split(':')
            rep_coord[rep_info[0]].append([int(line[0])+int(rep_info[1]), int(line[0])+int(rep_info[2])])
            rep_lens[rep_info[0]].append(int(rep_info[2])-int(rep_info[1]))
    print('Number of yL:' + str(len(rep_coord['yL'])))
    print('Number of yD:' + str(len(rep_coord['yD'])))
    print('Number of z2L:' + str(len(rep_coord['z2L'])))
    print('Number of z2D:' + str(len(rep_coord['z2D'])))
    #Plotting distribution of REPs length.
    fig=plt.figure(figsize=(15,15), dpi=100)
    #yL length
    plot0=plt.subplot2grid((2,2),(0,0), rowspan=1, colspan=1)     
    plot0.hist(rep_lens['yL'], color='#7FCE79', edgecolor='black', alpha=0.8)
    plot0.annotate('Mean yL length='+str(round(np.mean(rep_lens['yL']),2)), xy=(0.45, 0.9), xycoords='axes fraction', size=15)
    plot0.set_xlabel('yL length', size=17)
    plot0.set_ylabel('Number of yL REPs', size=17)
    plot0.set_title('yL', size=18)   
    #yD length
    plot1=plt.subplot2grid((2,2),(0,1), rowspan=1, colspan=1)     
    plot1.hist(rep_lens['yD'], color='#7FCE79', edgecolor='black', alpha=0.8)
    plot1.annotate('Mean yD length='+str(round(np.mean(rep_lens['yD']),2)), xy=(0.45, 0.9), xycoords='axes fraction', size=15)
    plot1.set_xlabel('yD length', size=17)
    plot1.set_ylabel('Number of yD REPs', size=17)
    plot1.set_title('yD', size=18)   
    #z2L length
    plot2=plt.subplot2grid((2,2),(1,0), rowspan=1, colspan=1)     
    plot2.hist(rep_lens['z2L'], color='#7FCE79', edgecolor='black', alpha=0.8)
    plot2.annotate('Mean z2L length='+str(round(np.mean(rep_lens['z2L']),2)), xy=(0.45, 0.9), xycoords='axes fraction', size=15)
    plot2.set_xlabel('z2L length', size=17)
    plot2.set_ylabel('Number of z2L REPs', size=17)
    plot2.set_title('z2L', size=18)   
    #z2D length
    plot3=plt.subplot2grid((2,2),(1,1), rowspan=1, colspan=1)     
    plot3.hist(rep_lens['z2D'], color='#7FCE79', edgecolor='black', alpha=0.8)
    plot3.annotate('Mean z2D length='+str(round(np.mean(rep_lens['z2D']),2)), xy=(0.45, 0.9), xycoords='axes fraction', size=15)
    plot3.set_xlabel('z2D length', size=17)
    plot3.set_ylabel('Number of z2D REPs', size=17)
    plot3.set_title('z2D', size=18)    
    plt.tight_layout()
    plt.savefig(outpath + "REPs_length_distributions.png", dpi=400, figsize=(15, 15)) 
    plt.close()        
    reps_in.close()
    return rep_coord

#######
#Returns REPs sequences using coordinates specified.
#######

def return_rep_seq(sequence, rep_coord, outpath):
    fileout=open(outpath+"REPs_sequences.fasta", 'w')
    fileout_in_so=open(outpath+"REPs_sequences_same_orientation.fasta", 'w')
    for rep_type, rep_set in rep_coord.items():
        for rep_c in rep_set:
            fileout.write('>'+rep_type+'_'+str(rep_c[0])+'_'+str(rep_c[1])+'\n')
            fileout.write(sequence[rep_c[0]-1:rep_c[1]]+'\n')            
            if rep_type in ['yD', 'z2L']:
                fileout_in_so.write('>'+rep_type+'_'+str(rep_c[1])+'_'+str(rep_c[0])+'_rc'+'\n')
                fileout_in_so.write(str(Seq(sequence[rep_c[0]-1:rep_c[1]]).reverse_complement())+'\n')
            else:
                fileout_in_so.write('>'+rep_type+'_'+str(rep_c[0])+'_'+str(rep_c[1])+'\n')
                fileout_in_so.write(sequence[rep_c[0]-1:rep_c[1]]+'\n')                
    fileout.close()
    fileout_in_so.close()   
    return

#######
#Returns "consensus" score for the set of intervals.
#######

def consensus_score_track(sp_rep_len, sp_rep_score):
    max_len=max(sp_rep_len)  
    consensus_score=[]
    for i in range(max_len): #Itarates positions
        pos_score=[]
        for j in range(len(sp_rep_score)): #Iterates intervals
            if i<len(sp_rep_score[j]):
                pos_score.append(sp_rep_score[j][i])
            else:
                continue
        consensus_score.append(np.mean(pos_score)) #Calculates mean score for current position j
    return consensus_score

#######
#Returns REPs scores using coordinates specified.
#######

def return_rep_score(genome_score, rep_coord, outpath):
    #Dict of lists which contain score traks correspondes to REP intervals.
    rep_score={'y': [], 'z2': []}
    #Dict of lists which contain lengths of REP intervals.
    rep_len={'y': [], 'z2': []}
    print('yL: ' + str(len(rep_coord['yL'])))
    print('yD: ' + str(len(rep_coord['yD'])))
    print('z2L: ' + str(len(rep_coord['z2L'])))
    print('z2D: ' + str(len(rep_coord['z2D'])))
    for rep_type, rep_set in rep_coord.items(): #Iterates REP types
        for rep_c in rep_set: #Itarates REP intervals
            print(rep_type)
            if rep_type=='yL':
                rep_score['y'].append(genome_score[rep_c[0]-1:rep_c[1]-1])
                rep_len['y'].append(len(rep_score['y'][-1]))
            elif rep_type=='yD':
                int_of_int=genome_score[rep_c[0]-1:rep_c[1]-1]
                int_of_int_inv=int_of_int[::-1]
                rep_score['y'].append(int_of_int_inv)
                rep_len['y'].append(len(rep_score['y'][-1]))
            elif rep_type=='z2L':
                int_of_int=genome_score[rep_c[0]-1:rep_c[1]-1]
                int_of_int_inv=int_of_int[::-1]                
                rep_score['z2'].append(int_of_int_inv)
                rep_len['z2'].append(len(rep_score['z2'][-1]))
            elif rep_type=='z2D':
                rep_score['z2'].append(genome_score[rep_c[0]-1:rep_c[1]-1])  
                rep_len['z2'].append(len(rep_score['z2'][-1]))
    print('Mean y REP len: ' + str(np.mean(rep_len['y'])))
    print('Mean z2 REP len: ' + str(np.mean(rep_len['z2'])))
    #"Consensus" score traks - contain mean score for every position of REP "alignment".
    y_score_cons=consensus_score_track(rep_len['y'], rep_score['y']) 
    z2_score_cons=consensus_score_track(rep_len['z2'], rep_score['z2'])
    rep_consensus_score={'y': y_score_cons[:34], 'z2': z2_score_cons[:36]}
    
    #Plotting REPs consensus score.
    fig=plt.figure(figsize=(15,15), dpi=100)    
    bin_width=1
    x_grid=np.arange(5, 36, 5)
    #Consensus for z2
    x_pos=np.arange(1, 37, 1)
    plot0=plt.subplot2grid((1,2),(0,0), rowspan=1, colspan=1) 
    plot0.set_xlim(1, 36)
    plot0.set_xticks(x_grid, minor=True)
    plot0.xaxis.grid(True, which='minor', linewidth=1, linestyle='--', color='black')
    plot0.bar(x_pos, z2_score_cons[:36], bin_width, color='#7FCE79', linewidth=1, edgecolor='black', align='edge')
    plot0.set_xlabel('z2 REP position', size=17)
    plot0.set_ylabel('Mean score', size=17)
    plot0.set_title('z2 REP', size=18)  
    #Consensus for y
    x_pos=np.arange(1, 35, 1)
    plot1=plt.subplot2grid((1,2),(0,1), rowspan=1, colspan=1) 
    plot1.set_xlim(1, 35)
    plot1.set_xticks(x_grid, minor=True)
    plot1.xaxis.grid(True, which='minor', linewidth=1, linestyle='--', color='black')    
    plot1.bar(x_pos, y_score_cons[:34], bin_width, color='#7FCE79', linewidth=1, edgecolor='black', align='edge')
    plot1.set_xlabel('y REP position', size=17)
    plot1.set_ylabel('Mean score', size=17)
    plot1.set_title('y REP', size=18) 
    plt.tight_layout()
    plt.savefig(outpath + "REPs_consensus_score.png", dpi=400, figsize=(15, 15)) 
    plt.close()      
    return rep_consensus_score

#######
#GCSs localization within classified REPs.
#######

def locate_gcs(GCSs_sets_dict, rep_coord, outpath):
    GCSs_in_REPs={}
    for rep_type, reps in rep_coord.items(): #Iterate REP type
        GCSs_in_REPs[rep_type]={}
        if rep_type in ['yL', 'z2D']:
            for rep in reps: #Iterate REPs
                for a, s in GCSs_sets_dict.items(): #Iterate GCSs sets
                    if a not in GCSs_in_REPs[rep_type]:
                        GCSs_in_REPs[rep_type][a]=[]
                    for gcs, info in s.items(): #Iterate GCSs
                        if rep[1]>gcs>rep[0]-3:
                            GCSs_in_REPs[rep_type][a].append(gcs-rep[0]+1)
                            print(len(GCSs_in_REPs[rep_type][a]))
        elif rep_type in ['yD', 'z2L']:
            for rep in reps: #Iterate REPs
                for a, s in GCSs_sets_dict.items(): #Iterate GCSs sets
                    if a not in GCSs_in_REPs[rep_type]:
                        GCSs_in_REPs[rep_type][a]=[]
                    for gcs, info in s.items(): #Iterate GCSs
                        if rep[1]>gcs>rep[0]-3:
                            GCSs_in_REPs[rep_type][a].append(rep[1]-gcs-5+1)
                            print(len(GCSs_in_REPs[rep_type][a]))
    GCSs_in_REPs_comb={'y':{}, 'z2':{}}
    for ab, gcs_set in GCSs_sets_dict.items():
        GCSs_in_REPs_comb['y'][ab]=sorted(GCSs_in_REPs['yL'][ab]+GCSs_in_REPs['yD'][ab])
        print(ab+' y REPs\' GCSs: ' + str(len(GCSs_in_REPs_comb['y'][ab])))
        print(ab+' yL REPs\' GCSs: ' + str(len(GCSs_in_REPs['yL'][ab])))
        print(ab+' yD REPs\' GCSs: ' + str(len(GCSs_in_REPs['yD'][ab])))
        GCSs_in_REPs_comb['z2'][ab]=sorted(GCSs_in_REPs['z2L'][ab]+GCSs_in_REPs['z2D'][ab])
        print(ab+' z2 REPs\' GCSs: ' + str(len(GCSs_in_REPs_comb['z2'][ab])))
        print(ab+' z2L REPs\' GCSs: ' + str(len(GCSs_in_REPs['z2L'][ab])))
        print(ab+' z2D REPs\' GCSs: ' + str(len(GCSs_in_REPs['z2D'][ab])))
    
    #Plotting distribution of GCSs within REPs.
    fig=plt.figure(figsize=(15,15), dpi=100)
    bins_y=np.linspace(1, 34, 34)
    #Cfx GCSs within y REPs
    plot0=plt.subplot2grid((4,2),(0,1), rowspan=1, colspan=1)     
    plot0.hist(GCSs_in_REPs_comb['y']['Cfx'], bins_y, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot0.set_xlabel('y REP position', size=17)
    plot0.set_ylabel('Number of Cfx GCSs', size=17)
    plot0.set_title('y REP', size=18)  
    #RifCfx GCSs within y REPs
    plot1=plt.subplot2grid((4,2),(1,1), rowspan=1, colspan=1)     
    plot1.hist(GCSs_in_REPs_comb['y']['RifCfx'], bins_y, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot1.set_xlabel('y REP position', size=17)
    plot1.set_ylabel('Number of RifCfx GCSs', size=17)
    plot1.set_title('y REP', size=18)   
    #Micro GCSs within y REPs
    plot2=plt.subplot2grid((4,2),(2,1), rowspan=1, colspan=1)     
    plot2.hist(GCSs_in_REPs_comb['y']['Micro'], bins_y, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot2.set_xlabel('y REP position', size=17)
    plot2.set_ylabel('Number of Micro GCSs', size=17)
    plot2.set_title('y REP', size=18)   
    #Oxo GCSs within y REPs
    plot3=plt.subplot2grid((4,2),(3,1), rowspan=1, colspan=1)     
    plot3.hist(GCSs_in_REPs_comb['y']['Oxo'], bins_y, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot3.set_xlabel('y REP position', size=17)
    plot3.set_ylabel('Number of Oxo GCSs', size=17)
    plot3.set_title('y REP', size=18)       
    
    bins_z2=np.linspace(1, 36, 36)
    #Cfx GCSs within z2 REPs
    plot4=plt.subplot2grid((4,2),(0,0), rowspan=1, colspan=1)     
    plot4.hist(GCSs_in_REPs_comb['z2']['Cfx'], bins_z2, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot4.set_xlabel('z2 REP position', size=17)
    plot4.set_ylabel('Number of Cfx GCSs', size=17)
    plot4.set_title('z2 REP', size=18) 
    #RifCfx GCSs within z2 REPs
    plot5=plt.subplot2grid((4,2),(1,0), rowspan=1, colspan=1)     
    plot5.hist(GCSs_in_REPs_comb['z2']['RifCfx'], bins_z2, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot5.set_xlabel('z2 REP position', size=17)
    plot5.set_ylabel('Number of RifCfx GCSs', size=17)
    plot5.set_title('z2 REP', size=18) 
    #Micro GCSs within z2 REPs
    plot6=plt.subplot2grid((4,2),(2,0), rowspan=1, colspan=1)     
    plot6.hist(GCSs_in_REPs_comb['z2']['Micro'], bins_z2, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot6.set_xlabel('z2 REP position', size=17)
    plot6.set_ylabel('Number of Micro GCSs', size=17)
    plot6.set_title('z2 REP', size=18) 
    #Oxo GCSs within z2 REPs
    plot7=plt.subplot2grid((4,2),(3,0), rowspan=1, colspan=1)     
    plot7.hist(GCSs_in_REPs_comb['z2']['Oxo'], bins_z2, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot7.set_xlabel('z2 REP position', size=17)
    plot7.set_ylabel('Number of Oxo GCSs', size=17)
    plot7.set_title('z2 REP', size=18)     
    plt.tight_layout()
    plt.savefig(outpath + "GCSs_relative_positions_distributions.png", dpi=400, figsize=(15, 15)) 
    plt.close()           
    return GCSs_in_REPs_comb

#######
#Write localization-classification data.
#######

def write_loc_class_GCSs(GCSs_in_REPs_comb, outpath):
    fileout=open(outpath+"GCSs_coordinates_within_REPs.txt", "w")
    fileout.write('REP type\tCondition\tGCSs coordinates\n')
    for rep_type, cond in GCSs_in_REPs_comb.items():
        for ab, gcss in cond.items():
            fileout.write(rep_type + '\t' + ab)
            for gcs in gcss:
                fileout.write('\t' + str(gcs))
            fileout.write('\n')
    fileout.close()
    return

#######
#Functions wrapper.
#######

def wraps_functions(input_dict, input_reps_path, score_inpath, outpath, seq_path):
    genome_seq=obtain_seq(seq_path)
    GCSs_sets_dict=trusted_GCSs_parsing(input_dict)
    genome_score=score_data_parser(score_inpath, 'score')
    rep_coord=REPs_parsing(input_reps_path, outpath)
    rep_consensus_score=return_rep_score(genome_score, rep_coord, outpath)
    return_rep_seq(genome_seq[0], rep_coord, outpath)
    GCSs_in_REPs_comb=locate_gcs(GCSs_sets_dict, rep_coord, outpath)
    write_loc_class_GCSs(GCSs_in_REPs_comb, outpath)
    return

wraps_functions(path_to_GCSs_files, REPs_coord_input, Score_data_path, Outputpath, Genome_seq_path)

print('Script ended its work succesfully!') 