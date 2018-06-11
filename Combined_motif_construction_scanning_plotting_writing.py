###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes sets of trusted GCSs as input, filters GCSs with highest N3E,
#makes a combined set consists of these GCSs, returns sequences under them and constructs
#PSSM matrix by the way getting rid of antibiotic-specific bias at positions forming the
#cleavage site. Than the script scans a sequence of interest with the PSSM, 
#returns the results of scanning, plots combined motif and writes it in a GC% degenerate and 
#in a non-degenerate forms.
###############################################

#######
#Packages to be imported.
#######

import scipy
from random import shuffle
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet  import IUPAC 
import matplotlib.pyplot as plt

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Input data - GCSs, TAB.
path_to_GCSs_files={'Cfx': '',
                    'Micro': '',
                    'Oxo': ''}

#Path to the E. coli genome (source of sequences for PFM/PWM construction), FASTA.
Genome_seq_path=''

#Path to the sequence to be scanned, FASTA.
Target_seq_path=''

#Prefix of the output path.
Output_data_prefix=''

###############################################
#Motif construction and scanning sequences of interest with it.
###############################################

#######
#FASTA sequences parsing.
#######

def obtain_seq(seq_path):
    seq_oi=open(seq_path, 'r')
    for record in SeqIO.parse(seq_oi, "fasta"):
        sequence=str(record.seq)
    seq_oi.close()      
    return sequence

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
                GCSs_dict[int(line[0])]=float(line[1])
            else:
                continue
        filein.close()
        GCSs_sets_dict[k]=GCSs_dict
        print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(GCSs_dict)))
    return GCSs_sets_dict

#######
#Selects GCSs with top height, makes combined set of GCSs.
#######

def sorting_combining(GCSs_dict):
    threshold=522
    GCSs_set_for_motif={}
    for k, v in GCSs_dict.items():
        list_of_keys_sbv=sorted(v, key=v.get)
        list_of_keys_sbv_slice=list_of_keys_sbv[-threshold:]
        for i in list_of_keys_sbv_slice:
            if i not in GCSs_set_for_motif:
                GCSs_set_for_motif[i]=v[i]
    print("Number of Cfx Micro Oxo GCSs are using for combined motif construction: " + str(len(GCSs_set_for_motif)))  
    return GCSs_set_for_motif 

#######
#Extracts sequences using GCSs coordinates (GCSs_ar) from DNA seq (genomefam), computes PWM and scans sequence of interest (genomefas).
#Returns file with scores.
#######

def motif_construction_and_analysis(GCSs_dict, genomefam, genomefas, outpath):
    #Frozen variables.
    background={'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}
    win_width_l=63
    win_width_r=67
    #Returns sequences under the GCSs.
    seqs=[]
    for k, v in GCSs_dict.items():
        seq=genomefam[(int(k) - win_width_l + 1):(int(k) + win_width_r + 1)]
        seqs.append(seq) 
    print('Number of sequences are using for motif construction: ' + str(len(seqs)))
    #Substitutes nucleotides at antibiotic-biased central positions with letters 
    #to be expected by chance taking into account the background nucleotides frequencies.
    hard_copy=['A']*(int(background['A']*len(seqs))+3) + ['C']*(int(background['C']*len(seqs))) +  ['G']*(int(background['G']*len(seqs))) + ['T']*(int(background['T']*len(seqs)))
    print(len(hard_copy))
    for k in range(65, 69):
        soft_copy=[]
        soft_copy.extend(hard_copy)
        shuffle(soft_copy)
        l=0
        for j in range(0, len(seqs)):
            ls=list(seqs[j])
            ls[k]=soft_copy[l]
            l=l+1
            seqs[j]="".join(ls)    
    #Create PWM and PSSM.
    instances=[]
    for seq in seqs:
        instances.append(Seq(seq))
    m=motifs.create(instances)
    pwm=m.counts.normalize()
    pssm=pwm.log_odds(background)
    #Scans forward sequence
    test_seq=Seq(str(genomefas), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq) 
    outfile=open(outpath + '_scan_forward_with_combined_motif.txt', 'w')
    for i in range(len(whole_genome_scores)):
        outfile.write(str(i+62) + '\t' + str(whole_genome_scores[i]) + '\t'+ str(test_seq[i+62-1]) + '\n')
    outfile.close()
    #Scans reverse complement sequence
    test_seq_rc=Seq(str(genomefas), IUPAC.unambiguous_dna).reverse_complement()
    whole_genome_scores_rc=pssm.calculate(test_seq_rc)   
    outfile_rc=open(outpath + '_scan_rc_with_combined_motif.txt', 'w')
    for i in range(len(whole_genome_scores_rc)):
        outfile_rc.write(str(i+63) + '\t' + str(whole_genome_scores_rc[-i]) + '\t'+ str(test_seq_rc[-(i+63)]) + '\n')
    outfile_rc.close()     
    return

#######
#Wraps functions for motif construction and for scanning sequences of interest with it.
#######

def Wrapper_motif_construct_scan(Source_genome_path, Target_genome_path, GCSs_files_paths, outpath):
    Source_sequence=obtain_seq(Source_genome_path)
    Target_sequence=obtain_seq(Target_genome_path)
    GCSs_sets=trusted_GCSs_parsing(GCSs_files_paths)
    GCSs_for_motif=sorting_combining(GCSs_sets)
    motif_construction_and_analysis(GCSs_for_motif, Source_sequence, Target_sequence, outpath)
    return

Wrapper_motif_construct_scan(Genome_seq_path, Target_seq_path, path_to_GCSs_files, Output_data_prefix)

###############################################
#The motif plotting and writing.
###############################################

#######
#Returns list of DNA seqs under the GCSs. Seqs have len=win_width.
#Writes sequences under the GCSs to file.
#######

def return_seqs(GCS_coords_dict, win_range, genomefa, outpath): 
    fileout=open(outpath + '_combined_motif_sequences.fasta', 'w')
    seqs=[]
    for k, v in GCS_coords_dict:
        seq=genomefa[int(int(k) - win_range[0] + 1):int(int(k) + win_range[1] + 1)]
        seqs.append(seq)
        fileout.write('>'+str(k)+'\n'+str(seq)+'\n')
    fileout.close()
    print('Number of sequences (GCSs) to analyse: ' + str(len(seqs)))
    return seqs

#######
#PFM construction.
#Scans sequences stack by columns, counts the number of particular letters.
#Returns a range of PFMs - "positional frequencies matrices" .
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
#Edits PPM in the cleavage site - substitutes real frequencies with background values.
#######

def gc_matrix_edit(matrix, background_freq):
    background_gc=float(background_freq['C']) + float(background_freq['G'])
    matrix[84]=background_gc
    matrix[85]=background_gc
    matrix[86]=background_gc
    matrix[87]=background_gc
    return matrix

#######
#Plots edited motif.
#######

def Plotting(matrix, matrix_type, win_width, outpath):
    x_axis=[]
    for i in range(len(matrix)):
        x_axis.append(-(win_width/2)+1+i)
    ax_range=[-win_width/2, win_width/2, 0.3, 0.95]
    plt.figure(figsize=(16, 6), dpi=900)
    plot1=plt.subplot()
    plot1.plot(x_axis, matrix, color='#7FCE79', linewidth=4, alpha=0.6)
    plot1.plot(x_axis, matrix, color='#454F24', linewidth=1, alpha=0.6)
    plot1.plot(x_axis, matrix, 'o', fillstyle='none', color='#7FCE79', markeredgecolor='#454F24', markersize=2, alpha=0.6)               
    plot1.tick_params(axis='both', direction='in', bottom='on', top='on', left='on', right='on')
    plot1.axis(ax_range)
    plot1.set_xlim(-win_width/2, win_width/2)
    plot1.set_xticks(np.concatenate((np.arange(-(win_width/2)+5, (win_width/2)+2, 10), [0, 3, -63, -17, 20, 66])))
    plot1.set_xlabel('Position, nt', size=17)
    plot1.set_ylabel(str(matrix_type), size=17)
    plt.show()
    plt.savefig(outpath + '_combined_motif_plot.png', dpi=600, figsize=(16, 6)) 
    return

#######
#Writes degenerate GC% PFM lacking antibiotic-biased frequencies for further Fourier analysis.
#######

def write_motif_deg(gc_matrix, win_width, outpath):
    fileout=open(outpath + '_deg_combined_motif_fourier.txt', 'w')
    fileout.write("#X\tY\n")
    for i in range(len(gc_matrix)):
        fileout.write(str((-win_width/2)+1+i) + '\t' + str(gc_matrix[i])+'\n')
    fileout.close()
    return

#######
#Writes non degenerate PFM edited.
#######

def write_line(ar, fileout):
    for i in range(len(ar)):
        fileout.write(str(ar[i])+'\t')
    fileout.write('\n')
    return

def write_motif(A, T, G, C, background, outpath):
    fileout=open(outpath + '_combined_motif_matrix.txt', 'w')
    A[84]=background['A']
    A[85]=background['A']
    A[86]=background['A']
    A[87]=background['A']  
    T[84]=background['T']
    T[85]=background['T']
    T[86]=background['T']
    T[87]=background['T']
    G[84]=background['G']
    G[85]=background['G']
    G[86]=background['G']
    G[87]=background['G']
    C[84]=background['C']
    C[85]=background['C']
    C[86]=background['C']
    C[87]=background['C']  
    write_line(A, fileout)
    write_line(T, fileout)
    write_line(G, fileout)
    write_line(C, fileout)
    fileout.close()
    return

#######
#Wraps functions for motif editing, plotting and writing.
#######

def Wrapper_motif_plotting_write(GCSs_dict, Source_genome_path, outpath):
    background={'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}
    window_width=170
    win_range=[(window_width/2)-2, (window_width/2)+2]
    Source_sequence=obtain_seq(Source_genome_path)
    GCSs_seqs=return_seqs(GCSs_dict, win_range, Source_sequence, outpath)
    pfms=make_PFM(GCSs_seqs)
    matrix_deg_red=gc_matrix_edit(pfms['GC'], background)
    Plotting(matrix_deg_red, 'GC', window_width, outpath)
    write_motif_deg(matrix_deg_red, window_width, outpath)
    write_motif(pfms['A'], pfms['T'], pfms['G'], pfms['C'], background, outpath)
    return

Wrapper_motif_plotting_write(path_to_GCSs_files, Genome_seq_path, Output_data_prefix)

print('Script ended its work succesfully!') 