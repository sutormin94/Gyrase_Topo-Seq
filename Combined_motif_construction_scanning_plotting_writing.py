###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes sets of trusted GCSs as input, filters GCSs with highest N3E,
#makes a combined set consists of these GCSs, returns sequences under them and constructs
#PSSM matrix by the way getting rid of antibiotic-specific bias at positions forming the
#cleavage site. Than the script scans a sequence of interest with the PSSM, 
#returns the results of scanning (WIG file), plots combined motif and writes it in a GC% degenerate and 
#in a non-degenerate forms.
###############################################

#######
#Packages to be imported.
#######

import os
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

#Input: GCSs data, TAB.
path_to_GCSs_files={'Cfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\Cfx_10mkM_trusted_GCSs.txt",
                    'Micro': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\Micro_trusted_GCSs.txt",
                    'Oxo': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\Oxo_trusted_GCSs.txt"}

#Input: path to the E. coli genome (source of sequences for PFM/PWM construction), FASTA.
Genome_seq_path="C:\Sutor\science\DNA-gyrase\Genomes\E_coli_w3110_G_Mu.fasta"

#Input: path to the sequence to be scanned, FASTA.
Target_seq_path="C:\Sutor\science\DNA-gyrase\Genomes\E_coli_w3110_G_Mu.fasta"
#Input: name of the target sequence ready to be scanned.
Target_seq_name="E_coli_w3110_G_Mu"
#Input: dataset name for WIG header.
Dataset_name="Cfx_Micro_Oxo_1828"

#Output: prefix of the output path.
Output_data_prefix="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\Combined_motif\\"
if not os.path.exists(Output_data_prefix):
    os.makedirs(Output_data_prefix)
    
#Output: path to the output WIG file.
path_to_res_score_files="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\Score_tracks\\"
if not os.path.exists(path_to_res_score_files):
    os.makedirs(path_to_res_score_files)
Output_score_wig=path_to_res_score_files + Target_seq_name + "_score.wig"


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
    #Estimate the threshold for the number of GCSs to take from each GCSs set to construct combined motif. 
    GCSs_sets_len=[]
    for k, v in GCSs_dict.items():
        GCSs_sets_len.append(len(v))
    threshold=min(GCSs_sets_len)
    #Return GCSs for combined motif construction.
    GCSs_set_for_motif={}
    for k, v in GCSs_dict.items():
        list_of_keys_sbv=sorted(v, key=v.get)
        list_of_keys_sbv_slice=list_of_keys_sbv[-threshold:]
        for i in list_of_keys_sbv_slice:
            if i not in GCSs_set_for_motif:
                GCSs_set_for_motif[i]=v[i]
    print("Number of Cfx Micro Oxo GCSs for combined motif construction: " + str(len(GCSs_set_for_motif)))  
    return GCSs_set_for_motif 

#######
#Extracts sequences using GCSs coordinates (GCSs_ar) from DNA seq (genomefam), computes PWM and scans sequence of interest (genomefas).
#Returns file with scores.
#######

def motif_construction_and_analysis(GCSs_dict, genomefam, genomefas, target_name, outpath):
    #Frozen variables.
    background={'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}
    win_width_l=63
    win_width_r=67
    #Returns sequences under the GCSs.
    seqs=[]
    for k, v in GCSs_dict.items():
        seq=genomefam[(int(k)-1-win_width_l):(int(k)-1+win_width_r)]
        seqs.append(seq) 
    print('Number of sequences for motif construction: ' + str(len(seqs)))
    print('Len of sequences for motif construction: ' + str(len(seqs[0])))
    #Substitutes nucleotides at antibiotic-biased central positions with letters 
    #to be expected by chance taking into account the background nucleotides frequencies.
    hard_copy=['A']*(int(background['A']*len(seqs))+3) + ['C']*(int(background['C']*len(seqs))) +  ['G']*(int(background['G']*len(seqs))) + ['T']*(int(background['T']*len(seqs)))
    print(len(hard_copy))
    for k in range(65, 69):
        soft_copy=[]
        soft_copy.extend(hard_copy)
        shuffle(soft_copy)
        for j in range(len(seqs)):
            ls=list(seqs[j])
            ls[k]=soft_copy[j]
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
    outfile=open(outpath + target_name + '_scan_forward_with_combined_motif.txt', 'w')
    dict_f={}
    #If test_seq len is equal to pssm, pssm.calculate() returns float32 but not the list.
    if len(test_seq)==win_width_l+win_width_r:
        one_pos_score=[]
        one_pos_score.append(whole_genome_scores)
        whole_genome_scores=one_pos_score
        dict_f[63+2]=[whole_genome_scores, test_seq[63+1]]
    for i in range(len(whole_genome_scores)):
        outfile.write(str(i+63+2) + '\t' + str(whole_genome_scores[i]) + '\t'+ str(test_seq[i+63+1]) + '\n')
        dict_f[i+63+2]=[whole_genome_scores[i], test_seq[i+63+1]]
    outfile.close()
    #Scans reverse complement sequence
    test_seq_rc=Seq(str(genomefas), IUPAC.unambiguous_dna).reverse_complement()
    whole_genome_scores_rc=pssm.calculate(test_seq_rc)   
    outfile_rc=open(outpath + target_name + '_scan_rc_with_combined_motif.txt', 'w')
    dict_rc={}
    #If test_seq_rc len is equal to pssm, pssm.calculate() returns float32 but not the list.
    if len(test_seq_rc)==win_width_l+win_width_r:
        one_pos_score_rc=[]
        one_pos_score_rc.append(whole_genome_scores_rc)
        whole_genome_scores_rc=one_pos_score_rc 
        dict_rc[67-1]=[whole_genome_scores_rc, test_seq_rc[67-1]]
    for i in range(len(whole_genome_scores_rc)):
        outfile_rc.write(str(i+67-1) + '\t' + str(whole_genome_scores_rc[-i-1]) + '\t'+ str(test_seq_rc[-(i+67-1)]) + '\n')
        dict_rc[i+67-1]=[whole_genome_scores_rc[-i-1], test_seq_rc[-(i+67-1)]]
    outfile_rc.close()     
    return dict_f, dict_rc


#######
#Combines fw and rc results of scanning procedure, prepares score WIG file.
#######

def combine_score_fw_rc(dict_f, dict_rc, genomefas, path_out, ds_name, chr_name):
    print("Length of the fw scan track: " + str(len(dict_f)))
    print("Length of the rc scan track: " + str(len(dict_rc)))
    #Preparing output WIG file.
    fileout=open(path_out, "w")
    fileout.write("track type=wiggle_0 name=\""+str(ds_name)+"\" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom=" + str(chr_name) + " start=1 step=1\n")       
    #Make combined score track.
    ar_max=[]
    for i in range(len(genomefas)):
        if i+1 in dict_f and i+1 in dict_rc:
            ar_max.append(max(dict_f[i+1][0], dict_rc[i+1][0]))
        elif i+1 in dict_f and i+1 not in dict_rc:
            ar_max.append(dict_f[i+1][0])
        elif i+1 not in dict_f and i+1 in dict_rc:
            ar_max.append(dict_rc[i+1][0])
        elif i+1 not in dict_f and i+1 not in dict_rc:
            ar_max.append(0)
    #Write score containing WIG file.
    for i in range(len(ar_max)):
        fileout.write(str(ar_max[i]) + "\n")
    fileout.close()
    return ar_max

#######
#Colocalization of top-scored genome sites with GCSs.
#######

def top_scored_sites_and_GCSs(ar_max, GCSs_sets_dict, outpath):
    fileout_prec=open(outpath + 'GCSs_matches_top_scored_genome_sites_precisely.txt', 'w')
    fileout_vic=open(outpath + 'GCSs_fall_into_5bp_vicinity_of_top_scored_genome_sites.txt', 'w')
    #Makes list of tuples with coordinates from ar_max scores array.
    Score_ar=[]
    for i in range(len(ar_max)):
        Score_ar.append([i+1, ar_max[i]])
    print(i)
    #Sort list by score.
    print(Score_ar[0])
    Score_ar.sort(key=lambda x: x[1], reverse=True)
    print(Score_ar[0])
    #Merges GCSs dicts into one dict.
    all_GCSs_dict={}
    for set_type, gcs_set in GCSs_sets_dict.items():
        for GCSs in gcs_set:
            if GCSs not in all_GCSs_dict:
                all_GCSs_dict[GCSs]=[set_type]
            else:
                all_GCSs_dict[GCSs].append(set_type)
    #Finds GCSs nearby top-scored genome sites.
    fileout_prec.write('Position\tScore\tGCS\n')
    fileout_vic.write('Position\tScore\tGCSs\n')
    for top_score in Score_ar[:15]:
        fileout_prec.write(str(top_score[0]) + '\t' + str(top_score[1]))
        fileout_vic.write(str(top_score[0]) + '\t' + str(top_score[1]))
        for GCSs in all_GCSs_dict:
            if top_score[0]==GCSs+1:
                fileout_prec.write('\t' + str(GCSs) + ' ' + str(all_GCSs_dict[GCSs]))
            if GCSs+3+5>=top_score[0]>=GCSs-3-5:
                fileout_vic.write('\t' + str(GCSs) + ' ' + str(all_GCSs_dict[GCSs]))
        fileout_prec.write('\n')
        fileout_vic.write('\n')
    fileout_prec.close()    
    fileout_vic.close()
    return

#######
#Wraps functions for motif construction and for scanning sequences of interest with it.
#######

def Wrapper_motif_construct_scan(Source_genome_path, Target_genome_path, target_name, GCSs_files_paths, outpath, outpath_top_score, wig_path_out, ds_name):
    Source_sequence=obtain_seq(Source_genome_path)[0]
    Target_sequence, chr_name=obtain_seq(Target_genome_path)
    GCSs_sets=trusted_GCSs_parsing(GCSs_files_paths)
    GCSs_for_motif=sorting_combining(GCSs_sets)
    f_rc_scan=motif_construction_and_analysis(GCSs_for_motif, Source_sequence, Target_sequence, target_name, outpath)
    ar_max=combine_score_fw_rc(f_rc_scan[0], f_rc_scan[1], Target_sequence, wig_path_out, ds_name, chr_name)
    top_scored_sites_and_GCSs(ar_max, GCSs_sets, outpath_top_score)
    return GCSs_for_motif


Motif_defined_GSCs_dict=Wrapper_motif_construct_scan(Genome_seq_path, Target_seq_path, Target_seq_name, path_to_GCSs_files, Output_data_prefix, path_to_res_score_files, Output_score_wig, Dataset_name)

###############################################
#The motif plotting and writing.
###############################################

#######
#Returns list of DNA seqs under the GCSs. Seqs have len=win_width.
#Writes sequences under the GCSs to file.
#######

def return_seqs(GCS_coords_dict, win_range, genomefa, outpath): 
    fileout=open(outpath + 'Combined_motif_originated_sequences.fasta', 'w')
    seqs=[]
    for k, v in GCS_coords_dict.items():
        seq=genomefa[int(int(k) - win_range[0] - 1):int(int(k) + win_range[1] - 1)]
        seqs.append(seq)
        fileout.write('>'+str(k)+'\n'+str(seq)+'\n')
    fileout.close()
    print('Number of sequences (GCSs) to be analysed: ' + str(len(seqs)))
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

def gc_matrix_edit(matrix, background_freq, win_width):
    background_gc=float(background_freq['C']) + float(background_freq['G'])
    Seq_centre=int((win_width/2))
    matrix[Seq_centre-1]=background_gc
    matrix[Seq_centre]=background_gc
    matrix[Seq_centre+1]=background_gc
    matrix[Seq_centre+2]=background_gc
    return matrix

#######
#Plots edited motif.
#######

def Plotting(matrix, matrix_type, win_width, outpath):
    x_axis=[]
    for i in range(len(matrix)):
        x_axis.append(-(win_width/2)+1+i)
    ax_range=[-win_width/2, win_width/2, 0.3, 0.95]
    plt.figure(figsize=(16, 6), dpi=100)
    plot1=plt.subplot()
    plot1.set_xticks([0], minor=True)
    plot1.xaxis.grid(True, which='minor', linewidth=0.5, linestyle='--', alpha=1)     
    plot1.plot(x_axis, matrix, color='#7FCE79', linewidth=4, alpha=0.6)
    plot1.plot(x_axis, matrix, color='#454F24', linewidth=1, alpha=0.6)
    plot1.plot(x_axis, matrix, 'o', fillstyle='none', color='#7FCE79', markeredgecolor='#454F24', markersize=2, alpha=0.6)               
    plot1.tick_params(axis='both', direction='in', bottom='on', top='on', left='on', right='on')
    plot1.axis(ax_range)
    plot1.set_xlim(-win_width/2, win_width/2)
    plot1.set_xticks(np.concatenate((np.arange(-(win_width/2)+5, (win_width/2)+2, 10), [0, 3, -63, -17, 20, 66])))
    plot1.set_xlabel('Position, nt', size=17)
    plot1.set_ylabel(str(matrix_type + '%'), size=17)
    #plt.show()
    plt.savefig(outpath + 'Combined_motif_plot.png', dpi=400, figsize=(16, 6))
    plt.close()
    return

#######
#Writes degenerate GC% PFM lacking antibiotic-biased frequencies for further Fourier analysis.
#######

def write_motif_deg(gc_matrix, win_width, outpath):
    fileout=open(outpath + 'Degenerate_combined_motif_for_Fourier_analysis.txt', 'w')
    fileout.write("#X\tY\n")
    for i in range(len(gc_matrix)):
        fileout.write(str(int((-win_width/2)+1+i)) + '\t' + str(gc_matrix[i])+'\n')
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

def write_motif(A, T, G, C, background, win_width, outpath):
    #Writes PFM matrix.
    fileout=open(outpath + 'Combined_motif_PFM.txt', 'w')
    Seq_centre=int((win_width/2))
    A[Seq_centre-1]=background['A']
    A[Seq_centre]=background['A']
    A[Seq_centre+1]=background['A']
    A[Seq_centre+2]=background['A']  
    T[Seq_centre-1]=background['T']
    T[Seq_centre]=background['T']
    T[Seq_centre+1]=background['T']
    T[Seq_centre+2]=background['T']
    G[Seq_centre-1]=background['G']
    G[Seq_centre]=background['G']
    G[Seq_centre+1]=background['G']
    G[Seq_centre+2]=background['G']
    C[Seq_centre-1]=background['C']
    C[Seq_centre]=background['C']
    C[Seq_centre+1]=background['C']
    C[Seq_centre+2]=background['C']  
    write_line(A, fileout)
    write_line(T, fileout)
    write_line(G, fileout)
    write_line(C, fileout)
    fileout.close()
    #Writes consensus sequence.
    fileout_cs=open(outpath + 'Combined_motif_consensus_sequence.fasta', 'w')
    fileout_cs.write('>Combined_motif_consensus_sequence\n')
    alphabet=['A', 'T', 'G', 'C']
    for i in range(len(A)):
        posit_examined=[A[i], T[i], G[i], C[i]]
        index_consensus=max(range(len(posit_examined)), key=posit_examined.__getitem__)
        fileout_cs.write(alphabet[index_consensus])
    fileout_cs.close()   
    return

#######
#Wraps functions for motif editing, plotting and writing.
#######

def Wrapper_motif_plotting_write(GCSs_form_motif_dict, Source_genome_path, outpath):
    background={'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}
    window_width=170
    win_range=[(window_width/2)-2, (window_width/2)+2]
    Source_sequence=obtain_seq(Source_genome_path)[0]
    GCSs_seqs=return_seqs(GCSs_form_motif_dict, win_range, Source_sequence, outpath)
    pfms=make_PFM(GCSs_seqs)
    matrix_deg_red=gc_matrix_edit(pfms['GC'], background, window_width)
    Plotting(matrix_deg_red, 'GC', window_width, outpath)
    write_motif_deg(matrix_deg_red, window_width, outpath)
    write_motif(pfms['A'], pfms['T'], pfms['G'], pfms['C'], background, window_width, outpath)
    return

Wrapper_motif_plotting_write(Motif_defined_GSCs_dict, Genome_seq_path, Output_data_prefix)

print('Script ended its work succesfully!') 