###############################################
##Dmitry Sutormin, 2024##
##Topo-Seq analysis. Calling of cleavage sites (CSs) strand-specifically. Milti-chromosome referense genomes are supported.##

#The script takes tetrades of WIG files contain N3E or N5E values: A+IP+, A+IP-, A-IP+, A-IP-.
#It smooths A+IP- and A-IP- tracks and divides A+IP+ and A-IP+ by them.
#Once obtains A+IP+_div and A-IP+_div the script performs Audic-Clavery
#statistic test and returns regions of A+IP+_div where i and i+5 positions are
#significantly higher than corresponding in A-IP+. These regions till now are called GCSs.
#GCSs are stored in the output TXT file. 
#Also two plots are generated: 1) signal coverage over the genome for treated and untreated samples;
#2) Motif expected to be under the GCSs.

#Requirements: TAB file with deletions.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count
import matplotlib.patheffects as PathEffects

#######
#Variables to be defined.
#######

#Path to the working directory
pwd="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\"
#Path to the file with regions to be omitted (e.g. deletions).
Deletions="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\Additional_genome_features\Deletions.broadPeak"
#Paths to the WIG files contain N3E or N5E that forms a tetrade: A+IP+, A+IP-, A-IP+, A-IP-.
#Name of the set (e.g., Cfx, RifCfx, and so on).

#A. thaliana gyrase Cfx data.
Tetrade_1={'A+IP+_N3E_F': pwd + "\WIG\\Topo_Seq_1_S5_edt_N3E_F.wig", 
           'A+IP+_N3E_R': pwd + "\WIG\\Topo_Seq_1_S5_edt_N3E_R.wig", 
           'A+IP-_N3E_F': pwd + "\WIG\\Topo_Seq_Mock_1_S7_edt_N3E_F.wig",
           'A+IP-_N3E_R': pwd + "\WIG\\Topo_Seq_Mock_1_S7_edt_N3E_R.wig",
           'A-IP+_N3E_F': pwd + "\WIG\\Topo_Seq_Mock_1_S7_edt_N3E_F.wig", 
           'A-IP+_N3E_R': pwd + "\WIG\\Topo_Seq_Mock_1_S7_edt_N3E_R.wig",
           'A-IP-_N3E_F': pwd + "\WIG\\Topo_Seq_Mock_1_S7_edt_N3E_F.wig",
           'A-IP-_N3E_R': pwd + "\WIG\\Topo_Seq_Mock_1_S7_edt_N3E_R.wig",
           'A+IP+_for':   pwd + "\WIG\\Topo_Seq_1_S5_edt_forward_depth.wig",
           'A+IP+_rev':   pwd + "\WIG\\Topo_Seq_1_S5_edt_reverse_depth.wig",
           'Tetrade name': 'Gyrase_Cfx_1'
           }

Tetrade_2={'A+IP+_N3E_F': pwd + "\WIG\\Topo_Seq_2_S6_edt_N3E_F.wig", 
           'A+IP+_N3E_R': pwd + "\WIG\\Topo_Seq_2_S6_edt_N3E_R.wig", 
           'A+IP-_N3E_F': pwd + "\WIG\\Topo_Seq_Mock_2_S8_edt_N3E_F.wig",
           'A+IP-_N3E_R': pwd + "\WIG\\Topo_Seq_Mock_2_S8_edt_N3E_R.wig",
           'A-IP+_N3E_F': pwd + "\WIG\\Topo_Seq_Mock_2_S8_edt_N3E_F.wig", 
           'A-IP+_N3E_R': pwd + "\WIG\\Topo_Seq_Mock_2_S8_edt_N3E_R.wig",
           'A-IP-_N3E_F': pwd + "\WIG\\Topo_Seq_Mock_2_S8_edt_N3E_F.wig",
           'A-IP-_N3E_R': pwd + "\WIG\\Topo_Seq_Mock_2_S8_edt_N3E_R.wig",
           'A+IP+_for':   pwd + "\WIG\\Topo_Seq_2_S6_edt_forward_depth.wig",
           'A+IP+_rev':   pwd + "\WIG\\Topo_Seq_2_S6_edt_reverse_depth.wig",
           'Tetrade name': 'Gyrase_Cfx_2'
           }


#Path to the reference genome
Genome="C:\\Users\sutor\OneDrive\\ThinkPad_working\Sutor\Science\Arabidopsis_gyrase\Sequencing_results\\GCA_000001735.2_TAIR10.1_genomic.fna\\GCA_000001735.2_TAIR10.1_genomic.fna"
#CSs calling p-value threshold.
P_value_threshold='0_05'
#Output folder
Path_for_output_1=os.path.join(pwd, 'GCS_calling', f'{Tetrade_1["Tetrade name"]}_GCSs_calling_{P_value_threshold}') 
Path_for_output_2=os.path.join(pwd, 'GCS_calling', f'{Tetrade_2["Tetrade name"]}_GCSs_calling_{P_value_threshold}') 

def check_path(Path_for_output):
    if not os.path.exists(Path_for_output):
        os.makedirs(Path_for_output)
    return 

check_path(Path_for_output_1)
check_path(Path_for_output_2)


#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    
    #Reads FASTA file with the reference genome. 
    chr_dict={}
    
    for record in SeqIO.parse(genome_path, "fasta"):
        chr_id=record.id
        chr_seq=str(record.seq)
        
        chr_dict[chr_id]=chr_seq
        print(f'Identified chromosome {chr_id}, length {len(chr_seq)} nt')
            
    return chr_dict

#######
#Opens and reads BED file with deletions coordinates.
#Example:
#GenomeID\tStart\tEnd
#NC_007779.1_w3110_Mu\t274500\t372148
#######

def deletions_info(del_path):
    
    del_dict={}
    filein=open(del_path, 'r')
    
    for line in filein:
        line=line.rstrip().split('\t')
        
        chr_id=line[0]
        
        if chr_id not in del_dict:
            del_dict[chr_id]=[int(line[1]), int(line[2])]
            
        else:
            del_dict[chr_id].append([int(line[1]), int(line[2])])
        
    filein.close()
    
    return del_dict

#######
#Returns nearby NE value if current position falls into deleted region of the genome.
#######

def get_value(i, ends, deletions):
    if i<0: #coordinate is out of the left genome border (start)
        j=len(ends)+i
    elif i>=len(ends): #coordinate is out of the right genome border (end)
        j=i-len(ends)
    else: #coordinate is within the genome borders
        check_in_del=0
        for dl in deletions: #check if coordinate falls into deletion
            if dl[1]>=i>=dl[0]:
                j=dl[1]-dl[0]+i+1
                check_in_del=1
        if check_in_del==0:
            j=i
    return ends[j]

#######
#Parses WIG file with N3/5E values.
#Computes a total number of Ends.
#######

def wig_parsing(wigfile):
    
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    
    NE_values_dict={}

    NE_values=[]
    Total_NE=0
    
    for line in wigin:
        line=line.rstrip().split(' ')
        
        if line[0]!='track':
        
            if line[0]=='fixedStep':
                
                if len(NE_values)>0:
                    
                    NE_values_dict[chr_id]=[NE_values, Total_NE]
                    
                    print(f'Total number of ends for chromosome {chr_id}: {Total_NE}')
                
                chr_id=line[1].split('=')[1]
                NE_values=[]
                Total_NE=0                
                
            else:
                
                NE_values.append(int(line[0]))
                Total_NE+=int(line[0])
            
    NE_values_dict[chr_id]=[NE_values, Total_NE]
    
    print(f'Total number of ends for chromosome {chr_id}: {Total_NE}')
    
    wigin.close()
    
    return NE_values_dict

#######
#Returns smoothed N3/5E tracks.
#Smoothing using sliding window (default - 200000 nt).
#######

def Smoothing(ends, deletions):
    smoothed=[]
    #Calculating the value for the first genome position
    mean=0.0
    window=100000
    window_float=float(window)
    for i in range(-window, window):
        mean=mean + get_value(i, ends, deletions)
    mean=mean/(2*window_float)
    smoothed.append(mean)
    #Calculating values for the part of the genome remains
    for i in range(1, len(ends)):
        mean=mean + (get_value(i+window, ends, deletions) - get_value(i-window, ends, deletions))/(2*window_float)
        smoothed.append(mean)
    return smoothed

#######
#Returns A+IP+/smoothed(A+IP-) and A-IP+/smoothed(A-IP-) tracks ready for GCSs calling.
#######

def norm_smooth_devide(treated_experiment_F, treated_experiment_R, treated_control_F, treated_control_R, 
                       untreated_experiment_F, untreated_experiment_R, untreated_control_F, untreated_control_R, 
                       deletions):
    
    #Normalization on the reads number
    #Adds pseudocounts to avoid zero values
    Min_total_NE=min(treated_experiment_F[1]+treated_experiment_R[1], treated_control_F[1]+treated_control_R[1], 
                     untreated_experiment_F[1]+untreated_experiment_R[1], untreated_control_F[1]+untreated_control_R[1])
    print('Min_total_NE: ' + str(Min_total_NE))
    
    treated_experiment_F_norm=[1.0 * (x) * Min_total_NE/treated_experiment_F[1] for x in treated_experiment_F[0]] #+A+IP norm F
    treated_control_F_norm=[1.0 * (x) * Min_total_NE/treated_control_F[1] for x in treated_control_F[0]] #+A-IP norm F
    untreated_experiment_F_norm=[1.0 * (x) * Min_total_NE/untreated_experiment_F[1] for x in untreated_experiment_F[0]] #-A+IP norm F
    untreated_control_F_norm=[1.0 * (x) * Min_total_NE/untreated_control_F[1] for x in untreated_control_F[0]] #-A-IP norm F
    treated_experiment_R_norm=[1.0 * (x) * Min_total_NE/treated_experiment_R[1] for x in treated_experiment_R[0]] #+A+IP norm R
    treated_control_R_norm=[1.0 * (x) * Min_total_NE/treated_control_R[1] for x in treated_control_R[0]] #+A-IP norm R
    untreated_experiment_R_norm=[1.0 * (x) * Min_total_NE/untreated_experiment_R[1] for x in untreated_experiment_R[0]] #-A+IP norm R
    untreated_control_R_norm=[1.0 * (x) * Min_total_NE/untreated_control_R[1] for x in untreated_control_R[0]] #-A-IP norm R
    
    #Control samples smoothing: A+IP- and A-IP-
    un_experiment_F_norm_sm=Smoothing(untreated_experiment_F_norm, deletions) #-A+IP F norm sm 
    un_control_F_norm_sm=Smoothing(untreated_control_F_norm, deletions) #-A-IP F norm sm
    un_experiment_R_norm_sm=Smoothing(untreated_experiment_R_norm, deletions) #-A+IP R norm sm 
    un_control_R_norm_sm=Smoothing(untreated_control_R_norm, deletions) #-A-IP R norm sm    
    
    #Pairwise division: +A+IP/-A+IP and +A-IP/-A-IP
    ends_divide_IP_F=[] #+A+IP F/-A+IP F
    ends_divide_IP_R=[] #+A+IP R/-A+IP R
    ends_divide_mock_F=[] #+A-IP F/-A-IP F
    ends_divide_mock_R=[] #+A-IP R/-A-IP R
    
    for i in range (len(treated_experiment_F_norm)):
        if un_experiment_F_norm_sm[i]!=0:
            ends_divide_IP_F.append(treated_experiment_F_norm[i]/un_experiment_F_norm_sm[i])
        else:
            ends_divide_IP_F.append((treated_experiment_F_norm[i]+1)/(un_experiment_F_norm_sm[i]+1))
        
        if un_experiment_R_norm_sm[i]!=0:
            ends_divide_IP_R.append(treated_experiment_R_norm[i]/un_experiment_R_norm_sm[i])
        else:
            ends_divide_IP_R.append((treated_experiment_R_norm[i]+1)/(un_experiment_R_norm_sm[i]+1))   
            
        if un_control_F_norm_sm[i]!=0:
            ends_divide_mock_F.append(treated_control_F_norm[i]/un_control_F_norm_sm[i])
        else:
            ends_divide_mock_F.append((treated_control_F_norm[i]+1)/(un_control_F_norm_sm[i]+1)) 
            
        if un_control_R_norm_sm[i]!=0:
            ends_divide_mock_R.append(treated_control_R_norm[i]/un_control_R_norm_sm[i])
        else:
            ends_divide_mock_R.append((treated_control_R_norm[i]+1)/(un_control_R_norm_sm[i]+1))
            
    Output_tracks=[ends_divide_IP_F, ends_divide_IP_R, ends_divide_mock_F, ends_divide_mock_R, 
                   un_experiment_F_norm_sm, un_experiment_R_norm_sm, un_control_F_norm_sm, un_control_R_norm_sm]
            
    return Output_tracks

#######
#Audic & Claverie statistics: borders of the confidential intervals (p-value=0.05, two-tailed test).
#From Audic & Claverie, 1997
#######

def AC_stat(x, P_value_threshold):
    x+=0
    #Confidential intervals borders (from Audic & Claverie, 1997).
    confidence=P_value_threshold
    if confidence=='0_05':
        AU_test=[5,7,9,11,12,14,16,17,19,20,22,23,24,26,27,28,30,31,32,34,35]
        AU_test20=20*1.75
        AU_test25=25*1.64
        AU_test30=30*1.60
        AU_test40=40*1.50
        AU_test50=50*1.44
        AU_test75=75*1.36   
        AU_test100=100*1.30
    elif confidence=='0_01':
        AU_test=[7,10,12,14,16,18,19,21,23,24,26,27,29,30,32,33,35,36,38,39,40]
        AU_test20=20*2
        AU_test25=25*1.88
        AU_test30=30*1.80
        AU_test40=40*1.68
        AU_test50=50*1.60
        AU_test75=75*1.48  
        AU_test100=100*1.40     
    #Estimation of a confidential interval higher border according to the value given - x.
    #Returns the interval border.
    if x<len(AU_test):
        int_border=AU_test[int(x)]
    elif 25>x>=20:
        int_border=AU_test20
    elif 30>x>=25:
        int_border=AU_test25
    elif 40>x>=30:
        int_border=AU_test30
    elif 50>x>=40:
        int_border=AU_test40
    elif 75>x>=50:
        int_border=AU_test50
    elif 100>x>=75:
        int_border=AU_test75
    else:
        int_border=AU_test100
    return int_border

#######
#Mark CSs in the genome.
#######

def mark_CSs(marked_css, dsCSs_ar, adsCSs_F_ar, adsCSs_R_ar, masking_type, gap_len):
    
    if masking_type=='precise':
        mask_len_l=0
        mask_len_r=gap_len+2
    elif masking_type=='extended':
        mask_len_l=2
        mask_len_r=gap_len+2+2
    
    for CS in dsCSs_ar:
        marked_css[CS-mask_len_l:CS+mask_len_r]+=np.array([1]*(mask_len_l+mask_len_r))
        
    for CS in adsCSs_F_ar:
        marked_css[CS-mask_len_l:CS+mask_len_r]+=np.array([1]*(mask_len_l+mask_len_r))
        
    for CS in adsCSs_R_ar:
        marked_css[CS-mask_len_l:CS+mask_len_r]+=np.array([1]*(mask_len_l+mask_len_r))

    return marked_css

#######
#Clean list of CSs; exclude overlapping CSs.
#######

def clean_CSs_ar(Marked_SCs, CSs_ar):
    
    CSs_ar_overlapped=[]
    CSs_ar_clean=[]
    
    for CS in CSs_ar:
        if max(Marked_SCs[CS:CS+5+1])>1:
            CSs_ar_overlapped.append(CS)
        else:
            CSs_ar_clean.append(CS)

    return CSs_ar_clean, CSs_ar_overlapped

#######
#Clean list of CSs; exclude CSs with satellite under-the-threshold CSs.
#######

def clean_sat_CSs(CSs_ar, IP_norm_div_F, IP_norm_div_R, mock_norm_div_F, mock_norm_div_R, P_value_threshold):
    
    CSs_ar_clean=[] 
    CSs_ar_sat=[]
    
    for CS in CSs_ar:
        if IP_norm_div_F[CS+5+2]<AC_stat(mock_norm_div_F[CS+5+2], P_value_threshold) and IP_norm_div_R[CS+5+2]<AC_stat(mock_norm_div_R[CS+5+2], P_value_threshold): 
            CSs_ar_clean.append(CS)
        else:
            CSs_ar_sat.append(CS)
               
    return CSs_ar_clean, CSs_ar_sat

#######
#Collect sequences under the CSs for motif plotting.
#######

def get_seq_for_motif(CS_ar, genome_fasta, win_width, orient, chr_id):
    seqs=[]
    if len(CS_ar)>0:
        for i in range(len(CS_ar)):
            if orient=='no':
                seq=genome_fasta[int(CS_ar[i])-int(win_width/2):int(CS_ar[i])+int(win_width/2)]
                seqs.append(seq)
            elif orient=='F':
                seq=genome_fasta[int(CS_ar[i])-int(win_width/2):int(CS_ar[i])+int(win_width/2)]
                seqs.append(seq)
            elif orient=='R':
                seq=genome_fasta[int(CS_ar[i])-int(win_width/2):int(CS_ar[i])+int(win_width/2)]
                seqs.append(seq)
            
    print(f'Number of sequences obtained for chromosome {chr_id}: {len(seqs)}')
    
    if len(seqs)==0:
        print(f'No CSs were called for chromosome {chr_id}!')
        
    return seqs

#######
#Plots the motif.
#######

def plot_the_motif(fname, seqs, genome_fasta, path_out):
    #PWM construction
    #Scans sequences stack by columns, counts the number of particular letters
    #Returns PFM (positional frequency matrix) - "matrix"
    matrix=[]
    template=seqs[0]
    for i in range(len(template)):
        column=[0, 0, 0, 0] #Corresponds to ['A', 'T', 'G', 'C']
        for j in range(len(seqs)):
            if seqs[j][i] in ['A', 'a']:
                column[0]+=1
            elif seqs[j][i] in ['T', 't']:
                column[1]+=1
            elif seqs[j][i] in ['G', 'g']:
                column[2]+=1
            elif seqs[j][i] in ['C', 'c']:
                column[3]+=1
        matrix.append(column)

    #Counts different combinations of nucleotides, only GC is used subsequently by default.
    #Returns degenerate PFMs.
    GC_percent = []
    GA_percent = []
    GT_percent = []
    AT_percent = []
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
        AT = float((int(matrix[i][0]) + int(matrix[i][1]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GA = float((int(matrix[i][0]) + int(matrix[i][2]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        CT = float((int(matrix[i][1]) + int(matrix[i][3]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        A = float((int(matrix[i][0]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        T = float((int(matrix[i][1]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        G = float((int(matrix[i][2]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        C = float((int(matrix[i][3]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GC_percent.append(GC)
        GT_percent.append(GT)
        AT_percent.append(AT)
        CT_percent.append(CT)
        GA_percent.append(GA)
        A_percent.append(A)
        T_percent.append(T)
        G_percent.append(G)
        C_percent.append(C)

    #GC statistics module
    #Counts average GC% over the whole genome
    GC_genome=GC_count(genome_fasta)/100
    print('GC% of the reference genome: ' + str(GC_genome))

    #Counts GC% p-value in the particular pwm column.
    #Returns p-value array and auxiliary Zero array for plotting.
    alignment_thick=len(seqs)
    pValue=[]
    Zero=[]
    for i in range(len(GC_percent)):
        pValue.append(scipy.stats.binom_test(float(GC_percent[i]) * alignment_thick, n=alignment_thick, p=GC_genome))
        Zero.append(1)

    #Plotting
    x_axis=[]
    for i in range(len(GC_percent)):
        x_axis.append(-111 + i)
    ax_range=[-110, +110, 0.2, 1]
    plt.figure(figsize=(16, 8), dpi=100)
    #GC% pwm plotting
    plt.suptitle(fname, fontsize=20)
    plot1=plt.subplot()
    plot1.plot(x_axis, GC_percent, color='green', linewidth=1)
    plot1.plot(x_axis, GC_percent, 'go', markersize=3)
    plot1.set_xticks(np.arange(-120, 112, 10))
    plot1.axis(ax_range)
    plot1.set_xlim(-110, 110)
    plot1.set_xticks([0], minor=True)
    plot1.xaxis.grid(True, which='minor', linewidth=0.5, linestyle='--', alpha=1)    
    plot1.annotate('GC%', xytext=(80, 0.65), xy=(40, 0.85), color='green', weight="bold", size=15)
    txt=plot1.annotate('p-value', xytext=(80, 0.60), xy=(-105, 0.64), color='cyan', weight="bold", size=15)
    txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='black')])   
    plot1.set_xlabel('Position, nt', size=17)
    plot1.set_ylabel('GC%', size=17)
    #p-value plotting
    plot2=plot1.twinx()
    plot2.plot(x_axis, pValue, 'k', linewidth=0.5, alpha=0.6)
    plot2.fill_between(x_axis, pValue, Zero, color='cyan', alpha=0.2)
    plot2.set_yticks(np.arange(0, 1.01, 0.01), minor=False)
    plot2.set_yscale('log')
    plot2.set_yticks([0.005], minor=True)
    plot2.yaxis.grid(True, which='minor', linewidth=1, linestyle='--', alpha=0.8)
    plot2.annotate('Confidence level = 0.005', xytext=(60, 0.0025), xy=(40, 0.8), color='black', size=15)
    plot2.set_ylim(0.0000001, 1.0)
    plot2.set_xlim(-110, 110)
    plot2.set_ylabel('p-value, logarithmic scale', size=17)
    plt.show()
    plt.savefig(path_out+fname+'_motif_expected.png', dpi=300, figsize=(16, 8))
    plt.close()
    return

#######
#Prepare the metasignal enrichment at CSs.
#######

def prepare_metasignal_data(CSs_ar, IP_norm_div, win_width):
    
    CSs_metasignal=np.array([0.0]*(win_width+win_width+5+1))
    for CS_coord in CSs_ar:
        if (CS_coord-win_width>0) and (CS_coord+win_width+5+1)<len(IP_norm_div):
            Signal_around_CS=np.array(IP_norm_div[CS_coord-win_width:CS_coord+win_width+5+1])
            CSs_metasignal+=Signal_around_CS
        
    CSs_metasignal=CSs_metasignal/len(CSs_ar)    
    
    return CSs_metasignal

#######
#Prepare the combined metasignal enrichment at CSs for F and R CSs.
#######

def prepare_metasignal_data_FR(CSs_ar_F, CSs_ar_R, IP_norm_div, win_width):
    
    CSs_metasignal=np.array([0.0]*(win_width+win_width+5+1))
    for CS_coord in CSs_ar_F:
        if (CS_coord-win_width>0) and (CS_coord+win_width+5+1)<len(IP_norm_div):
            Signal_around_CS=np.array(IP_norm_div[CS_coord-win_width:CS_coord+win_width+5+1])
            CSs_metasignal+=Signal_around_CS
            
    for CS_coord in CSs_ar_R:
        if (CS_coord-win_width>0) and (CS_coord+win_width+5+1)<len(IP_norm_div):
            Signal_around_CS=np.array(list(np.array(IP_norm_div[CS_coord-win_width:CS_coord+win_width+5+1]))[::-1])
            CSs_metasignal+=Signal_around_CS
        
    CSs_metasignal=CSs_metasignal/(len(CSs_ar_F)+len(CSs_ar_R)) 
    
    return CSs_metasignal

#######
#Plots the metasignal enrichment at CSs.
#######

def plot_metasignal_CSs(fname, CSs_type, CSs_metasignal, IP_norm_div, plot_type, win_width, path_out):    
    
    positions=np.arange(-win_width, win_width+5+1, 1)
    plt.figure(figsize=(10, 6), dpi=100)
    plot1=plt.subplot(111)   
    if plot_type=='smoothed':
        plot1.plot(positions, CSs_metasignal,   linestyle='-',  color='#6e4bc2', linewidth=2.5, alpha=1, label=f'{fname} over {CSs_type}', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        plot1.fill_between(positions, [0]*len(CSs_metasignal), CSs_metasignal, facecolor='#6e4bc2', alpha=0.4, interpolate=True)  
        ticks=np.arange(-win_width,win_width+5+1, 20).tolist()   
        plot1.set_xticks(ticks, minor=False)
        ticks_lables=ticks
        ticks_lables[ticks.index(0)]=CSs_type   
        plot1.axvline(0, color='black', linestyle='--', alpha=0.7, linewidth=1)
    elif plot_type=='barplot':
        plot1.bar(positions, CSs_metasignal,   width=1,  align='center', color='#6e4bc2', alpha=1, label=f'{fname} over {CSs_type}', zorder=6) #Def linewidth=1.5; #R123 -0; R12 +0.03; R3 -0.07
        ticks_min=np.arange(-win_width,win_width+5+1, 1).tolist()   
        plot1.set_xticks(ticks_min, labels=[], minor=True) 
        ticks_coords=list(np.arange(-win_width, 1, 2))+list(np.arange(1, win_width+5+1, 2))
        ticks_labels=list(np.arange(-win_width-1, 0, 2))+list(np.arange(1, win_width+5+1, 2))         
        plot1.set_xticks(ticks_coords, labels=ticks_labels, minor=False)
        plot1.axvline(0.5, color='black', linestyle='--', alpha=0.7, linewidth=1)
        plot1.axvline(4.5, color='black', linestyle='--', alpha=0.7, linewidth=1)

    plot1.legend(fontsize=15, frameon=False)    
    plot1.set_xlabel('Distance, bp', size=20)
    plot1.set_ylabel(f'FE', size=20)   
    plot1.set_ylim([-0.05*max(CSs_metasignal), 1.1*max(CSs_metasignal)])
    if plot_type=='smoothed':
        plt.savefig(f'{path_out}\\{fname}_FE_over_{CSs_type}_{win_width}bp.png', dpi=400, figsize=(10, 6), transparent=False)  #Def size - 10, 6; Closer look - 3, 6
        plt.savefig(f'{path_out}\\{fname}_FE_over_{CSs_type}_{win_width}bp.svg', dpi=400, figsize=(10, 6), transparent=False)  #Def size - 10, 6; Closer look - 3, 6
    elif plot_type=='barplot':
        plt.savefig(f'{path_out}\\{fname}_N3E_over_{CSs_type}_{win_width}bp.png', dpi=400, figsize=(10, 6), transparent=False)  #Def size - 10, 6; Closer look - 3, 6
        plt.savefig(f'{path_out}\\{fname}_N3E_over_{CSs_type}_{win_width}bp.svg', dpi=400, figsize=(10, 6), transparent=False)  #Def size - 10, 6; Closer look - 3, 6
    plt.show()
    plt.close()     
    
    return

#######
#Plots the enrichment signal over the genome: +A+IP/smoothed(-A+IP) and +A-IP/smoothed(-A-IP)
#######

def plot_enrichment_signal(fname, IP_nd_ends, mock_nd_ends, un_IP_sm, un_mock_sm, deletions, path_out):
    #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
    #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    mpl.rcParams['agg.path.chunksize']=10000
    #Scaling smoothed tracks to make them visible on the plot.
    max_element=max(IP_nd_ends+mock_nd_ends) #Max N3E value of experimental tracks
    max_element_IP_sm=max(un_IP_sm)
    max_element_mock_sm=max(un_mock_sm)
    un_IP_sm=[(max_element/2)*x/max_element_IP_sm for x in un_IP_sm]
    un_mock_sm=[(max_element/2)*x/max_element_mock_sm for x in un_mock_sm]
    #Regions to be masked (e.g. deletions).  
    mask_array=[]
    for k in range(len(IP_nd_ends)):
        check_in_del=0
        for dl in deletions:
            if dl[1]>=k>=dl[0]:
                mask_array.append(True)
                check_in_del=1
        if check_in_del==0:
            mask_array.append(False)
    IPed=np.ma.masked_array(IP_nd_ends, mask=mask_array)
    mock=np.ma.masked_array(mock_nd_ends, mask=mask_array)
    un_IPed=np.ma.masked_array(un_IP_sm, mask=mask_array)
    un_mock=np.ma.masked_array(un_mock_sm, mask=mask_array)
    #Plotting the distribution of the signal around the genome for IPed and mock samples.
    xcoord=np.arange(0,len(IPed))
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle(fname, fontsize=20)
    plot1=plt.subplot() 
    plot1.plot(xcoord, IPed, '-', label='+A+IP/smoothed(-A+IP)', color='blue', linewidth=1)
    plot1.plot(xcoord, mock, '-', label='+A-IP/smoothed(-A-IP)', color='orange', linewidth=1)
    plot1.plot(xcoord, un_IPed, '-', label='smoothed(-A+IP)', color='#5bbdff', linewidth=3)
    plot1.plot(xcoord, un_mock, '-', label='smoothed(-A-IP)', color='#ed781f', linewidth=3)    
    plot1.set_xlabel('Genome position, nt', size=17)
    plot1.set_ylabel('Signal enrichment', size=17)
    plot1.legend(loc='upper right')
    plt.show()
    plt.savefig(path_out+fname+'_signal_enrichment.png', dpi=300, figsize=(16, 8))
    plt.close()
    return

#######
#Write CSs coordinates in simplified tab and BroadPeak formats.
#######

def write_CSs_coordinates(IP_norm_div, CSs_ar, CS_type, Tet_ID, chr_id, path_out):
    
    if CS_type in ['ssCS_F', 'ssCS_R', 'All_ssCS']:
        GCs_broadpeak_out=open(f'{path_out}{Tet_ID}_raw_{CS_type}_called.BroadPeak', 'w')
        CSs_out=open(f'{path_out}{Tet_ID}_raw_{CS_type}_called.txt', 'w')        
    else:   
        GCs_broadpeak_out=open(f'{path_out}{Tet_ID}_raw_{CS_type}_{chr_id}_called.BroadPeak', 'w')
        CSs_out=open(f'{path_out}{Tet_ID}_raw_{CS_type}_{chr_id}_called.txt', 'w')
    
    CSs_out.write(f'{CS_type}_coordinate\tN3E\n')
    i=0
    for CS in CSs_ar:
        CSs_out.write(f'{CS+1}\t{IP_norm_div[CS]}\t{IP_norm_div[CS+5]}\n')
        if CS_type=='dsCS':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+6}\t{CS_type}_{i}\t{np.mean([IP_norm_div[CS], IP_norm_div[CS+5]])}\t.\t-1\t-1\t-1\n')
        elif CS_type=='adsCS_F':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+6}\t{CS_type}_{i}\t{IP_norm_div[CS]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='adsCS_R':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+6}\t{CS_type}_{i}\t{IP_norm_div[CS+5]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='All_dsCS':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+6}\t{CS_type}_{i}\t{np.mean([IP_norm_div[CS], IP_norm_div[CS+5]])}\t.\t-1\t-1\t-1\n')
        elif CS_type=='ssCS_F':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+6}\t{CS_type}_{i}\t{IP_norm_div[CS]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='ssCS_R':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+6}\t{CS_type}_{i}\t{IP_norm_div[CS+5]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='All_ssCS':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+6}\t{CS_type}_{i}\t{np.max([IP_norm_div[CS], IP_norm_div[CS+5]])}\t.\t-1\t-1\t-1\n')        
        elif CS_type=='dsCS_5':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+7}\t{CS_type}_{i}\t{np.mean([IP_norm_div[CS], IP_norm_div[CS+6]])}\t.\t-1\t-1\t-1\n')
        elif CS_type=='adsCS_5_F':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+7}\t{CS_type}_{i}\t{IP_norm_div[CS]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='adsCS_5_R':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+7}\t{CS_type}_{i}\t{IP_norm_div[CS+6]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='All_dsCS_5':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+7}\t{CS_type}_{i}\t{np.mean([IP_norm_div[CS], IP_norm_div[CS+6]])}\t.\t-1\t-1\t-1\n')   
        elif CS_type=='dsCS_6':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+8}\t{CS_type}_{i}\t{np.mean([IP_norm_div[CS], IP_norm_div[CS+7]])}\t.\t-1\t-1\t-1\n')
        elif CS_type=='adsCS_6_F':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+8}\t{CS_type}_{i}\t{IP_norm_div[CS]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='adsCS_6_R':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+8}\t{CS_type}_{i}\t{IP_norm_div[CS+7]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='All_dsCS_6':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+8}\t{CS_type}_{i}\t{np.mean([IP_norm_div[CS], IP_norm_div[CS+7]])}\t.\t-1\t-1\t-1\n')      
        elif CS_type=='dsCS_3':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+5}\t{CS_type}_{i}\t{np.mean([IP_norm_div[CS], IP_norm_div[CS+4]])}\t.\t-1\t-1\t-1\n')
        elif CS_type=='adsCS_3_F':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+5}\t{CS_type}_{i}\t{IP_norm_div[CS]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='adsCS_3_R':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+5}\t{CS_type}_{i}\t{IP_norm_div[CS+4]}\t.\t-1\t-1\t-1\n')
        elif CS_type=='All_dsCS_3':
            GCs_broadpeak_out.write(f'{chr_id}\t{CS}\t{CS+5}\t{CS_type}_{i}\t{np.mean([IP_norm_div[CS], IP_norm_div[CS+4]])}\t.\t-1\t-1\t-1\n')

        i+=1
        
    CSs_out.close()
    GCs_broadpeak_out.close()    
    
    return




#######
#Wraps all the functions: data normalization, GCSs calling, 
#plotting, writing the data.
#######

def GCSs_caller(tetrade_dictionary, deletions_inpath, genome_path, P_value_threshold, path_out):
    
    #Parsing deletions
    deletions=deletions_info(deletions_inpath)
    
    #Reading referense genome.
    genome_fasta=read_genome(genome_path)
    
    #Define samples within the tetrade.
    Tet_ID=tetrade_dictionary['Tetrade name']
    print('Now we are working with: ' + str(Tet_ID))
    
    #WIG parsing, total NE counting (for further normalization on reads number)
    treated_experiment_F=wig_parsing(tetrade_dictionary['A+IP+_N3E_F'])   #+A+IP F N3E
    treated_experiment_R=wig_parsing(tetrade_dictionary['A+IP+_N3E_R'])   #+A+IP R N3E
    treated_control_F=wig_parsing(tetrade_dictionary['A+IP-_N3E_F'])      #+A-IP F N3E
    treated_control_R=wig_parsing(tetrade_dictionary['A+IP-_N3E_R'])      #+A-IP R N3E
    untreated_experiment_F=wig_parsing(tetrade_dictionary['A-IP+_N3E_F']) #-A+IP F N3E
    untreated_experiment_R=wig_parsing(tetrade_dictionary['A-IP+_N3E_R']) #-A+IP R N3E
    untreated_control_F=wig_parsing(tetrade_dictionary['A-IP-_N3E_F'])    #-A-IP F N3E
    untreated_control_R=wig_parsing(tetrade_dictionary['A-IP-_N3E_R'])    #-A-IP R N3E
    treated_exp_for_cd=wig_parsing(tetrade_dictionary['A+IP+_for'])       #+A+IP F coverage depth
    treated_exp_rev_cd=wig_parsing(tetrade_dictionary['A+IP+_rev'])       #+A+IP R coverage depth
    
    for chr_id, chr_seq in genome_fasta.items():
        
        print(f'GCSs calling for chromosome {chr_id}')
        
        if chr_id in deletions:
            deletions_ar=deletions[chr_id]
        else:
            deletions_ar=[]
    
        #Obtain pairwise divided tracks: A+IP+/A+IP- and A-IP+/A-IP-.
        ends_fitting=norm_smooth_devide(treated_experiment_F[chr_id], treated_experiment_R[chr_id], treated_control_F[chr_id], treated_control_R[chr_id], 
                                        untreated_experiment_F[chr_id], untreated_experiment_R[chr_id], untreated_control_F[chr_id], untreated_control_R[chr_id], 
                                        deletions_ar)
    
        IP_norm_div_F=ends_fitting[0]
        IP_norm_div_R=ends_fitting[1]
        mock_norm_div_F=ends_fitting[2]
        mock_norm_div_R=ends_fitting[3]
        un_IP_sm_F=ends_fitting[4]
        un_IP_sm_R=ends_fitting[5]
        un_mock_sm_F=ends_fitting[6] 
        un_mock_sm_R=ends_fitting[7]
    
        ####
        #### Convential 4-bp dsCSs detection and visualization procedure.
        #### 
        
        #dsCSs calling using the Audic, Claverie test (Audic & Claverie, 1998).
        #Returns arrays of dsCSs, that contains coordinates of the left wall of the gap. 
        dsCSs_ar=[]
    
        for i in range(len(IP_norm_div_F)-1-5):
            
            if IP_norm_div_F[i]>AC_stat(mock_norm_div_F[i], P_value_threshold) and IP_norm_div_R[i+5]>AC_stat(mock_norm_div_R[i+5], P_value_threshold):
                #Additional verification of dsCSs using the coverage depth tracks.
                Gap_mean=np.mean(treated_exp_for_cd[chr_id][0][i+1:i+5])+np.mean(treated_exp_rev_cd[chr_id][0][i+1:i+5])
                Left_mean_F=np.mean(treated_exp_for_cd[chr_id][0][i-3:i+1])
                Left_mean_R=np.mean(treated_exp_rev_cd[chr_id][0][i-3:i+1])
                Right_mean_F=np.mean(treated_exp_for_cd[chr_id][0][i+5:i+9])
                Right_mean_R=np.mean(treated_exp_rev_cd[chr_id][0][i+5:i+9])
                Left_ss_clv=sum(IP_norm_div_F[i-1:i+2])
                Right_ss_clv=sum(IP_norm_div_R[i+4:i+7])
                
                if (Gap_mean<(Left_mean_F+Left_mean_R)) and (Gap_mean<(Right_mean_F+Right_mean_R)):  #A 4-bp gap exists.
                    
                    dsCSs_ar.append(i)
     
                    
        print(f'\nNumber of dsCSs found for chromosome {chr_id}: {len(dsCSs_ar)}')
    
        #Clean the CSs calling data. Exclude overlapping CSs.
        #Marked_SCs_template=np.array([0]*len(IP_norm_div_F))
        #Marked_SCs=mark_CSs(Marked_SCs_template, dsCSs_ar, adsCSs_F_ar, adsCSs_R_ar, 'precise', 4)
        #dsCSs_ar, dsCSs_ar_ovrlp=clean_CSs_ar(Marked_SCs, dsCSs_ar)
        #print(f'Number of non-overlapping dsCSs found for chromosome {chr_id}: {len(dsCSs_ar)}')        
    
        #Plotting: possible motif and signal distribution over the genome.
        #win_width - width of area of interest under GCSs center. 
        win_width=220

        #Returns "seqs" array contains sequences under the GCSs within the win_width vicinity of the GCSs.
        ds_seqs=get_seq_for_motif(dsCSs_ar, chr_seq, win_width, 'no', chr_id)
        
        #Plotting the motif to be expected.
        plot_the_motif(f'{Tet_ID}_dsCSs_{chr_id}', ds_seqs, chr_seq, path_out)
    
        #Plotting the metasignal enrichment at CSs.
        IP_norm_div=list(np.array(IP_norm_div_F)+np.array(IP_norm_div_R))
        dsCSs_meta=prepare_metasignal_data(dsCSs_ar, IP_norm_div, 10)
        plot_metasignal_CSs(f'{Tet_ID}_N3E_dsCSs_{chr_id}', 'dsCSs', dsCSs_meta, IP_norm_div, 'barplot', 10, path_out)
    
        treated_exp_cd=list(np.array(treated_exp_for_cd[chr_id][0])+np.array(treated_exp_rev_cd[chr_id][0]))
        dsCSs_meta=prepare_metasignal_data(dsCSs_ar, treated_exp_cd, 1000)
        plot_metasignal_CSs(f'{Tet_ID}_cd_dsCSs_{chr_id}', 'dsCSs', dsCSs_meta, treated_exp_cd, 'smoothed', 1000, path_out)
 
        #Plotting the distribution of the signal around the genome for treated and untreated.
        mock_norm_div=list(np.array(mock_norm_div_F)+np.array(mock_norm_div_R))
        un_IP_sm=list(np.array(un_IP_sm_F)+np.array(un_IP_sm_R))
        un_mock_sm=list(np.array(un_mock_sm_F)+np.array(un_mock_sm_R))
        plot_enrichment_signal(f'{Tet_ID}_{chr_id}', IP_norm_div, mock_norm_div, un_IP_sm, un_mock_sm, deletions_ar, path_out)    

        #Writes GCSs data.
        #Converts coordinates to 1-start system.
        write_CSs_coordinates(IP_norm_div, dsCSs_ar, 'dsCS', Tet_ID, chr_id, path_out)
    
    
    return


GCSs_caller(Tetrade_1, Deletions, Genome, P_value_threshold, Path_for_output_1)
GCSs_caller(Tetrade_2, Deletions, Genome, P_value_threshold, Path_for_output_2)


print('Script ended its work succesfully!')
