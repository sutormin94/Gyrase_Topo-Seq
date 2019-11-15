###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

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
pwd="/data/Gyrase/Data_preparation"
#Path to the file with regions to be omitted (e.g. deletions).
Deletions="/data/Gyrase/Genomes_tracks/Deletions.bed"
#Paths to the WIG files contain N3E or N5E that forms a tetrade: A+IP+, A+IP-, A-IP+, A-IP-.
#Name of the set (e.g. Cfx, RifCfx, Micro, Oxo and so on).
'''
'A+IP-': pwd + "/RifCfx/WIG/RifCfx_IN_Mu_122mkM_10mkM_1_edt_N3E.wig",
         'A-IP+': pwd + "/Un/WIG/Un_IP_Mu_2_edt_N3E.wig", 
         'A-IP-': pwd + "/RifCfx/WIG/Rif_IN_Mu_122mkM_1_edt_N3E.wig",
         'Tetrade name': 'RifCfx_1'
         }
'''
Tetrade={'A+IP+': pwd + "/Un/WIG/Un_IP_Mu_2_edt_N3E.wig", 
         'A+IP-': pwd + "/Un/WIG/Un_IN_Mu_2_edt_N3E.wig",
         'A-IP+': pwd + "/Un/WIG/Un_IP_Mu_1_edt_N3E.wig", 
         'A-IP-': pwd + "/Un/WIG/Un_IN_Mu_1_edt_N3E.wig",
         'Tetrade name': 'Un_2_1'
         }
'''
Tetrade={'A+IP+': "/data/Gyrase/Stuff/Pipe/Fragments_ends/Data_second_mapping/Cfx_IP_10_mkM_2_ends.wig", 
         'A+IP-': "/data/Gyrase/Stuff/Pipe/Fragments_ends/Data_second_mapping/Cfx_IN_10_mkM_3_ends.wig",
         'A-IP+': "/data/Gyrase/Stuff/Pipe/Fragments_ends/Data_second_mapping/Un_IN_1_ends.wig", 
         'A-IP-': "/data/Gyrase/Stuff/Pipe/Fragments_ends/Data_second_mapping/Un_IN_2_ends.wig",
         'Tetrade name': 'Cfx_3_Un_old'
         }
'''
#Path to the reference genome
Genome="/data/Gyrase/Genomes_tracks/E_coli_w3110_G_Mu.fasta"
#Output folder
Path_for_output=pwd + "/RifCfx/GCSs_calling_0_05_Un/" + Tetrade['Tetrade name'] + "/"
if not os.path.exists(Path_for_output):
    os.makedirs(Path_for_output)


#######
#Opens and reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the reference genome.
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genomefa=str(record.seq)
    return genomefa

#######
#Opens and reads BED file with deletions coordinates.
#Example:
#GenomeID\tStart\tEnd
#NC_007779.1_w3110_Mu\t274500\t372148
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return del_ar

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
    NE_values=[]
    Total_NE=0
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(int(line[0]))
            Total_NE+=int(line[0])
    print('Total number of ends: ' + str(Total_NE))
    wigin.close()
    return NE_values, Total_NE

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

def norm_smooth_devide(ex_file_path, cont_file_path, un_ex_file_path, un_cont_file_path, deletions):
    #WIG parsing, total NE counting (for further normalization on reads number)
    treated_experiment=wig_parsing(ex_file_path) #+A+IP
    treated_control=wig_parsing(cont_file_path) #+A-IP
    untreated_experiment=wig_parsing(un_ex_file_path) #-A+IP
    untreated_control=wig_parsing(un_cont_file_path) #-A-IP
    #Normalization on the reads number
    #Adds pseudocounts to avoid zero values
    Min_total_NE=min(treated_experiment[1], treated_control[1], 
                     untreated_experiment[1], untreated_control[1])
    print('Min_total_NE: ' + str(Min_total_NE))
    treated_experiment_norm=[1.0 * (x + 1) * Min_total_NE/treated_experiment[1] for x in treated_experiment[0]] #+A+IP norm
    treated_control_norm=[1.0 * (x + 1) * Min_total_NE/treated_control[1] for x in treated_control[0]] #+A-IP norm
    untreated_experiment_norm=[1.0 * (x + 1) * Min_total_NE/untreated_experiment[1] for x in untreated_experiment[0]] #-A+IP norm
    untreated_control_norm=[1.0 * (x + 1) * Min_total_NE/untreated_control[1] for x in untreated_control[0]] #-A-IP norm
    #Control samples smoothing: A+IP- and A-IP-
    un_experiment_norm_sm=Smoothing(untreated_experiment_norm, deletions) #-A+IP norm sm 
    un_control_norm_sm=Smoothing(untreated_control_norm, deletions) #-A-IP norm sm
    #Pairwise division: +A+IP/-A+IP and +A-IP/-A-IP
    ends_divide_IP=[] #+A+IP/-A+IP
    ends_divide_mock=[] #+A-IP/-A-IP
    for i in range (len(treated_experiment_norm)):
        if treated_experiment_norm[i]!=0 and un_experiment_norm_sm[i]!=0:
            ends_divide_IP.append(treated_experiment_norm[i]/un_experiment_norm_sm[i])
        else:
            ends_divide_IP.append(0)
        if treated_control_norm[i]!=0 and un_control_norm_sm[i]!=0:
            ends_divide_mock.append(treated_control_norm[i]/un_control_norm_sm[i])
        else:
            ends_divide_mock.append(0) 
    return ends_divide_IP, ends_divide_mock, un_experiment_norm_sm, un_control_norm_sm

#######
#Audic & Claverie statistics: borders of the confidential intervals (p-value=0.05, two-tailed test).
#From Audic & Claverie, 1997
#######

def AC_stat(x):
    x+=-1
    #Confidential intervals borders (from Audic & Claverie, 1997).
    confidence=0.05
    if confidence==0.05:
        AU_test=[5,7,9,11,12,14,16,17,19,20,22,23,24,26,27,28,30,31,32,34,35]
        AU_test20=20*1.75
        AU_test25=25*1.64
        AU_test30=30*1.60
        AU_test40=40*1.50
        AU_test50=50*1.44
        AU_test75=75*1.36   
        AU_test100=100*1.30
    elif confidence==0.01:
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
    xcoord=np.arange(0,4647999)
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
#Wraps all the functions: data normalization, GCSs calling, 
#plotting, writing the data.
#######

def GCSs_caller(tetrade_dictionary, deletions_inpath, genome_path, path_out):
    #Parsing deletions
    deletions=deletions_info(deletions_inpath)
    
    #Define samples within the tetrade.
    Tet_ID=tetrade_dictionary['Tetrade name']
    print('Now we are working with: ' + str(Tet_ID))
    treated_control=tetrade_dictionary['A+IP-']
    treated_experiment=tetrade_dictionary['A+IP+']
    untreated_control=tetrade_dictionary['A-IP-']
    untreated_experiment=tetrade_dictionary['A-IP+']
    
    #Obtain pairwise divided tracks: A+IP+/A+IP- and A-IP+/A-IP-.
    ends_fitting=norm_smooth_devide(treated_experiment, treated_control, untreated_experiment, untreated_control, deletions)
    IP_norm_div=ends_fitting[0]
    mock_norm_div=ends_fitting[1]
    un_IP_sm=ends_fitting[2]
    un_mock_sm=ends_fitting[3]    
    
    #GCSs calling procedure.
    #GCSs calling using Audic, Claverie test (Audic & Claverie, 1998).
    #Returns array of GCSs, that contains coordinates of the left wall of the gap.      
    GCSs=[]
    for i in range(len(IP_norm_div)-1-5):
        if IP_norm_div[i]>AC_stat(mock_norm_div[i]) and IP_norm_div[i+5]>AC_stat(mock_norm_div[i+5]):
            GCSs.append(i)
    print('Number of GCSs just found: ' + str(len(GCSs)))
    thr=25000
    if (len(GCSs)>thr):
        print('Number of GCSs is extremely high! The threshold is ' + str(thr) + '.\nJust warning...') 
    
    #Plotting: possible motif and signal distribution over the genome.
    #win_width - width of area of interest under GCSs center. 
    win_width=220
    genome_fasta=read_genome(genome_path)
    #Returns "seqs" array contains sequences under the GCSs within the win_width vicinity of the GCSs.
    seqs=[]
    for i in range(len(GCSs)):
        seq=genome_fasta[int(GCSs[i])-win_width/2:int(GCSs[i])+win_width/2]
        seqs.append(seq)
    print('Number of sequences obtained: ' + str(len(seqs)))
    if len(seqs)==0:
        print('No GCSs were called!')
    
    #Plotting the motif to be expected.
    plot_the_motif(Tet_ID, seqs, genome_fasta, path_out)
    
    #Plotting the distribution of the signal around the genome for treated and untreated.
    plot_enrichment_signal(Tet_ID, IP_norm_div, mock_norm_div, un_IP_sm, un_mock_sm, deletions, path_out)    

    #Writes GCSs data.
    #Converts coordinates to 1-start system.
    GCSs_out=open(path_out+Tet_ID+'_raw_GCSs_called.txt', 'w')
    GCSs_out.write('GCSs_coordinate\tN3E\n')
    for GCS in GCSs:
        GCSs_out.write(str(GCS+1) + '\t' + str(min(IP_norm_div[GCS], IP_norm_div[GCS+5])) + '\n')  
    GCSs_out.close()
    return

GCSs_caller(Tetrade, Deletions, Genome, Path_for_output)

print('Script ended its work succesfully!')
