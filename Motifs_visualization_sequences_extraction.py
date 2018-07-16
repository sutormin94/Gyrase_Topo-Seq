###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes sets of trusted GCSs as input and plots motifs using the sequences under the GCSs.
#Also it writes sequences and motif to the files.
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import SeqUtils

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Input data - GCSs, TAB.
path_to_GCSs_files={'Cfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\Cfx_10mkM_trusted_GCSs.txt",
                    'RifCfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\RifCfx_trusted_GCSs.txt",
                    'Micro': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\Micro_trusted_GCSs.txt",
                    'Oxo': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\Oxo_trusted_GCSs.txt"}

#Path to the genome FASTA.
Genome_path="C:\Sutor\science\DNA-gyrase\Genomes\E_coli_w3110_G_Mu.fasta"

#Path for the output.
Output_path="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\Motif\\"
if not os.path.exists(Output_path):
        os.makedirs(Output_path)

#######
#Trusted GCSs data parsing.
#######

def trusted_GCSs_parsing(input_dict):
        GCSs_sets_dict={}
        for k, v in input_dict.items():
                ar=[]
                filein=open(v, 'r')
                for line in filein:
                        line=line.rstrip().split('\t')
                        if line[0] not in ['GCSs_coordinate']:
                                ar.append(int(line[0]))
                        else:
                                continue
                filein.close()
                GCSs_sets_dict[k]=ar
                print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(ar)))
        return GCSs_sets_dict

#######
#Genome sequence parsing.
#######

def genome_seq(genome_path):
        genome=open(genome_path, 'r')
        for record in SeqIO.parse(genome, "fasta"):
                genome_sequence=str(record.seq)
        genome.close()
        print('Whole genome average GC: ' + str(SeqUtils.GC(genome_sequence)))
        print('Whole genome length: ' + str(len(genome_sequence)))        
        return genome_sequence

#######
#Returns list of DNA seqs under the GCSs. Seqs have len=win_width.
#Writes sequences under the GCSs to file.
#######

def return_seqs(GCS_coords, win_range, genomefa, filepath_full_len, filepath_6bp_LOGO): 
        fileout=open(filepath_full_len, 'w')
        fileout_6bp_LOGO=open(filepath_6bp_LOGO, 'w')
        seqs=[]
        for i in range(len(GCS_coords)):
                seq=genomefa[int(GCS_coords[i]- win_range[0] - 1):int(GCS_coords[i]+ win_range[1] - 1)]
                seq_6bp_LOGO=genomefa[int(GCS_coords[i]-1):int(GCS_coords[i]-1+6)]
                seqs.append(seq)
                fileout.write('>'+str(GCS_coords[i])+'\n'+str(seq)+'\n')
                fileout_6bp_LOGO.write('>'+str(GCS_coords[i])+'\n'+str(seq_6bp_LOGO)+'\n')
        fileout.close()
        fileout_6bp_LOGO.close()
        print('Number of sequences (GCSs) analysing: ' + str(len(seqs)))
        return seqs

#######
#PFM construction.
#Scans sequences stack by columns, counts the number of particular letters.
#Returns a range of PFMs - "positional frequencies matrixes" .
#######

def make_PFM(seqs_list):
        matrix=[]
        template=seqs_list[0]
        for i in range(len(template)):
                column=[0, 0, 0, 0]
                for j in range(len(seqs_list)):
                        if seqs_list[j][i] == str('A'):
                                column[0] = column[0] + 1
                        elif seqs_list[j][i] == str('T'):
                                column[1] = column[1] + 1
                        elif seqs_list[j][i] == str('G'):
                                column[2] = column[2] + 1
                        elif seqs_list[j][i] == str('C'):
                                column[3] = column[3] + 1
                matrix.append(column)
        #Returns a range of PFMs.
        GC_percent = []
        GT_percent = []
        CT_percent = []
        A_percent = []
        T_percent = []
        G_percent = []
        C_percent = []
        for i in range(len(matrix)):
                GC = float((int(matrix[i][2]) + int(matrix[i][3]))) / (
                        int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                GT = float((int(matrix[i][1]) + int(matrix[i][2]))) / (
                        int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                CT = float((int(matrix[i][1]) + int(matrix[i][3]))) / (
                        int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                A = float((int(matrix[i][0]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                T = float((int(matrix[i][1]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                G = float((int(matrix[i][2]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                C = float((int(matrix[i][3]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
                GC_percent.append(GC)
                GT_percent.append(GT)
                CT_percent.append(CT)
                A_percent.append(A)
                T_percent.append(T)
                G_percent.append(G)
                C_percent.append(C)
        return {'Num_seqs': len(seqs_list), 'A': A_percent, 'T': T_percent, 'G': G_percent, 'C': C_percent, 'CT': CT_percent, 'GT': GT_percent, 'GC': GC_percent}

#######
#Writes PFM data to file.
#######

def write_motif(ar, filepath, coord_shift):
        fileout=open(filepath, 'w')
        fileout.write("#X\tY\n")
        for i in range(len(ar)):
                fileout.write(str((-coord_shift/2)+1+i) + '\t' + str(ar[i])+'\n')
        fileout.close()
        return

#######
#Plotting the motif.
#Matrix type - type of the PFM to plot.
#######

def Plotting(PFMs_set, title, matrix_type, write_out, win_width):
        x_axis=[]
        for i in range(len(PFMs_set['Cfx'])):
                x_axis.append(-(win_width/2)+1+i)      
        print('len(x_axis)=' + str(len(x_axis)))
        ax_range = [-win_width/2, win_width/2, 0.35, 0.9]
        plt.figure(dpi=100, figsize=(16, 6))
        plt.suptitle(str(title), fontsize=20)
        plot1 = plt.subplot()
        plot1.set_xticks([0], minor=True)
        plot1.xaxis.grid(True, which='minor', linewidth=0.5, linestyle='--', alpha=1)            
        #Cfx
        plot1.plot(x_axis, PFMs_set['Cfx'], color='#7FCE79', linewidth=4, alpha=0.6)
        plot1.plot(x_axis, PFMs_set['Cfx'], color='#454F24', linewidth=1, alpha=0.6)
        plot1.plot(x_axis, PFMs_set['Cfx'], 'o', fillstyle='none', color='#7FCE79', markeredgecolor='#454F24', markersize=2, alpha=0.6)        
        #Rif_Cfx
        #plot1.plot(x_axis, PFMs_set['RifCfx'], color='#BAE85C', linewidth=4, alpha=0.6)
        #plot1.plot(x_axis, PFMs_set['RifCfx'], color='#4D590D', linewidth=1, alpha=0.6)
        #plot1.plot(x_axis, PFMs_set['RifCfx'], 'o', fillstyle='none', color='#BAE85C', markeredgecolor='#4D590D', markersize=2, alpha=0.6)       
        #Micro
        plot1.plot(x_axis, PFMs_set['Micro'], color='#ff878b', linewidth=4, alpha=0.7)
        plot1.plot(x_axis, PFMs_set['Micro'], color='#7D212B', linewidth=1, alpha=1)
        plot1.plot(x_axis, PFMs_set['Micro'], 'o', fillstyle='none', color='#ff878b', markeredgecolor='#7D212B', markersize=2, alpha=1)
        #Oxo
        plot1.plot(x_axis, PFMs_set['Oxo'], color='#8991ff', linewidth=4)
        plot1.plot(x_axis, PFMs_set['Oxo'], color='#470A59', linewidth=1)
        plot1.plot(x_axis, PFMs_set['Oxo'], 'o', fillstyle='none', color='#8991ff', markeredgecolor='#470A59', markersize=2, alpha=0.8)   
        #Tracks annotation
        plot1.annotate('Ciprofloxacin', xytext=(-75, 0.8), xy=(40, 0.85), color='#7FCE79', weight="bold", size=15)
        #plot1.annotate('Rifampicin Ciprofloxacin', xytext=(-75, 0.8), xy=(40, 0.85), color='#BAE85C', weight="bold", size=15)
        plot1.annotate('Microcin B17', xytext=(-75, 0.75), xy=(40, 0.85), color='#ff878b', weight="bold", size=15)
        plot1.annotate('Oxolinic acid', xytext=(-75, 0.70), xy=(40, 0.85), color='#8991ff', weight="bold", size=15)        
        #Set axis parameters
        plot1.tick_params(axis='both', direction='in', bottom='on', top='on', left='on', right='on')
        plot1.axis(ax_range)
        plot1.set_xlim(-win_width/2, win_width/2)
        plot1.set_xticks(np.concatenate((np.arange(-(win_width/2)+5, (win_width/2)+2, 10), [0, 3, -63, -17, 20, 66])))
        plot1.set_xlabel('Position, nt', size=17)
        plot1.set_ylabel(str(matrix_type), size=17)
        #plt.show()
        plt.savefig(write_out, dpi=400, figsize=(16, 6)) 
        plt.close()
        return

#######
#Wraps all the functions together.
#######

def wrap_function(GCSs_input, genome_input_path, output_path):
        win_width=170
        win_range=[(win_width/2)-2, (win_width/2)+2]
        PFM_type='GC'
        plot_title='Gyrase motif with different inhibitors'
        GCSs_dict=trusted_GCSs_parsing(GCSs_input)
        genome_sequence=genome_seq(genome_input_path)
        dict_of_PFMs={}
        for k, v in GCSs_dict.items():
                sequences_list=return_seqs(v, win_range, genome_sequence, output_path+str(k)+'_sequences_under_GCSs_full.fasta', output_path+str(k)+'_sequences_under_GCSs_6bp_LOGO.fasta')
                PFMs=make_PFM(sequences_list)
                write_motif(PFMs[PFM_type], output_path+str(k)+'_GC_pfm.txt', win_width)
                dict_of_PFMs[k]=PFMs[PFM_type]
        Plotting(dict_of_PFMs, plot_title, PFM_type, output_path+'Gyrase_motif_trusted_GCSs_Cfx_Micro_Oxo'+str(PFM_type)+'.png', win_width)
        return

wrap_function(path_to_GCSs_files, Genome_path, Output_path)

print('Script ended its work succesfully!') 