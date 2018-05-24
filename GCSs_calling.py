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
###############################################

#######
#Packages to be imported
#######

import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils import GC as GC_count

#######
#Variables to be defined
#######

print('Variables to be defined:')

#Path to the file with regions to be omitted (e.g. deletions).
Deletions=''
#Paths to the WIG files contain N3E or N5E that forms a tetrade: A+IP+, A+IP-, A-IP+, A-IP-.
#Name of the set (e.g. Cfx, RifCfx, Micro, Oxo and so on).
Tetrade={'A+IP+': '', 
         'A+IP-': '',
         'A-IP+': '', 
         'A-IP-': '',
         'Tetrade name': ''
         }
#Path to the reference genome
Genome=''
#Output folder
Path_for_output=''

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
#Opens and reads TAB file with deletions coordinates.
#Example:
#Deletions\tStart\tEnd
#1\t100000\t200000
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ['Deletions']:
            del_ar.append([line[1], line[2]])
    filein.close()
    return del_ar

List_of_deletions=deletions_info(Deletions)

#######
#Returns nearby NE value if current position falls into deleted region of the genome.
#######

def get_value(i, ends, deletions):
    if i<0:
        j=len(ends)+i
    elif i>=len(ends):
        j=i-len(ends)
    check_in_del=0
    for dl in deletions:
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

def wig_parsing(wigin):
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
    mean=0.0
    window=100000
    window_float=window*1.0
    for i in range(-1*window, window):
        mean=mean + get_value(i, ends, deletions)
    mean=mean/(2*window_float)
    smoothed.append(mean)
    for i in range(1, len(ends)):
        mean=mean + (get_value(i+window, ends, deletions) - get_value(i-window, ends, deletions))/(2*window_float)
        smoothed.append(mean)
    return smoothed

#######
#Returns A+IP+/smoothed(A+IP-) and A-IP+/smoothed(A-IP-) tracks ready for GCSs calling
#######

def norm_smooth_devide(ex_file_path, cont_file_path, un_ex_file_path, un_cont_file_path, deletions):
    #WIG parsing, total NE counting (for further normalization on reads number)
    treated_experiment=wig_parsing(ex_file_path)
    treated_control=wig_parsing(cont_file_path)
    untreated_experiment=wig_parsing(un_ex_file_path)
    untreated_control=wig_parsing(un_cont_file_path)
    #Normalization on the reads number
    Min_total_NE=min(treated_experiment[1], treated_control[1], untreated_experiment[1], untreated_control[1])
    treated_experiment_norm=[1.0 * x * Min_total_NE/treated_experiment[1] for x in treated_experiment[0]]
    treated_control_norm=[1.0 * x * Min_total_NE/treated_control[1] for x in treated_control[0]]
    untreated_experiment_norm=[1.0 * x * Min_total_NE/untreated_experiment[1]  for x in untreated_experiment[0]]
    untreated_control_norm=[1.0 * x * Min_total_NE/untreated_control[1] for x in untreated_control[0]]    
    #Control samples smoothing: A+IP- and A-IP-
    control_norm_sm=Smoothing(treated_control_norm, deletions)
    un_control_norm_sm=Smoothing(untreated_control_norm, deletions)
    #Adds pseudocounts to avoid zero values
    exper_ends=[1.0 * x + 1.0  for x in treated_experiment_norm] #A+IP+
    control_ends=[1.0 * x + 1.0  for x in control_norm_sm] #A+IP-
    un_exper_ends=[1.0 * x + 1.0  for x in untreated_experiment_norm] #A-IP+
    un_control_ends=[1.0 * x + 1.0 for x in un_control_norm_sm] #A-IP-
    #Pairwise division: A+IP+/A+IP- and A-IP+/A-IP-
    ends_divide_treated=[] #A+IP+/A+IP-
    ends_divide_untreated=[] #A-IP+/A-IP-
    for i in range (len(exper_ends)):
        if exper_ends[i]!=0 and control_ends[i]!=0:
            ends_divide_treated.append(exper_ends[i]/control_ends[i])
        else:
            ends_divide_treated.append(0)
        if un_exper_ends[i]!=0 and un_control_ends[i]!=0:
            ends_divide_untreated.append(un_exper_ends[i]/un_control_ends[i])
        else:
            ends_divide_untreated.append(0) 
    return ends_divide_treated, ends_divide_untreated

#######
#Audic & Claverie statistics: borders of the confidential intervals (p-value=0.05, two-tailed test).
#From Audic & Claverie, 1997
#######

def AC_stat(x):
    #Confidential intervals borders (from Audic & Claverie, 1997).
    AU_test=[5,7,9,11,12,14,16,17,19,20,22,23,24,26,27,28,30,31,32,34,35]
    AU_test20=20*1.75
    AU_test25=25*1.64
    AU_test30=30*1.60
    AU_test40=40*1.50
    AU_test50=50*1.44
    AU_test75=75*1.36
    AU_test100=100*1.30
    #Estimation of a confidential interval higher border according to the value given - x.
    #Returns the interval border.
    if len(AU_test)>x:
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
        x_axis.append(-109 + i)
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
    plot1.annotate('GC%', xytext=(85, 0.68), xy=(40, 0.85), color='green', weight="bold", size=15)
    plot1.annotate('p-value', xytext=(85, 0.63), xy=(-105, 0.64), color='cyan', weight="bold", size=15)
    plot1.set_xlabel('Position, nt', size=17)
    plot1.set_ylabel('GC%', size=17)
    #p-value plotting
    plot2=plot1.twinx()
    plot2.plot(x_axis, pValue, 'k', linewidth=0.5, alpha=0.6)
    plot2.fill_between(x_axis, pValue, Zero, color='cyan', alpha=0.2)
    plot2.set_yticks(np.arange(0, 1.01, 0.01), minor=False)
    plot2.set_yscale('log')
    plot2.set_yticks([0.005], minor=True)
    plot2.yaxis.grid(True, which='minor', linewidth=2, linestyle='--', alpha=0.3)
    plot2.annotate('Confidence level = 0.005', xytext=(66, 0.0032), xy=(40, 0.8), color='black', size=15)
    plot2.set_ylim(0.0000001, 1.0)
    plot2.set_xlim(-110, 110)
    plot2.set_ylabel('p-value, logarithmic scale', size=17)
    plt.show()
    plt.savefig(path_out+fname+'_motif_expected.png', dpi=300, figsize=(16, 8))
    plt.close()
    return

#######
#Plots the enrichment signal over the genome: A+IP+/smoothed(A+IP-) and A-IP+/smoothed(A-IP-)
#######

def plot_enrichment_signal(fname, treated_nd_ends, untreated_nd_ends, deletions, path_out):	
    #Regions to be masked (e.g. deletions).  
    mask_array=[]
    for k in range(len(treated_nd_ends)):
        check_in_del=0
        for dl in deletions:
            if dl[1]>=k>=dl[0]:
                mask_array.append(True)
                check_in_del=1
        if check_in_del==0:
            mask_array.append(False)
    treated=np.ma.masked_array(treated_nd_ends, mask=mask_array)
    unteated=np.ma.masked_array(untreated_nd_ends, mask=mask_array)
    #Plotting the distribution of the signal around the genome for treated and untreated.
    xcoord=np.arange(0,4647999)
    plt.figure(figsize=(16, 8), dpi=100)
    plt.plot(xcoord, treated, '.', label='A+IP+/smoothed(A+IP-)', color='blue')
    plt.plot(xcoord, unteated, '.', label='A-IP+/smoothed(A-IP-)', color='orange')
    plt.set_xlabel('Genome position, nt', size=17)
    plt.set_ylabel('Signal enrichment', size=17)
    plt.title(fname)
    plt.show()
    plt.savefig(path_out+fname+'_signal_enrichment.png', dpi=300, figsize=(16, 8))
    plt.close()
    return

#######
#Wraps all the functions: data normalization, GCSs calling, 
#plotting, writing the data.
#######

def GCSs_caller(tetrade_dictionary, deletions, genome_path, path_out):
    #Define samples within the tetrade.
    Tet_ID=tetrade_dictionary['Tetrade name']
    print('Now we are working with: ' + str(Tet_ID))
    treated_control=tetrade_dictionary['A+IP-']
    treated_experiment=tetrade_dictionary['A+IP+']
    untreated_control=tetrade_dictionary['A-IP-']
    untreated_experiment=tetrade_dictionary['A-IP+']
    
    #Obtain pairwise divided tracks: A+IP+/A+IP- and A-IP+/A-IP-.
    ends_fitting=norm_smooth_devide(treated_experiment, treated_control, untreated_experiment, untreated_control, deletions)
    treated_norm_div=ends_fitting[0]
    unterated_norm_div=ends_fitting[1]
    
    #GCSs calling procedure.
    #GCSs calling using Audic, Claverie test (Audic & Claverie, 1998).
    #Returns array of GCSs, that contains coordinates of the first position of a gap.
    #(or the position of the gap left wall in the case of 1-start-point coordinates).      
    GCSs=[]
    for i in range(len(treated_norm_div)-1-5):
        if treated_norm_div[i]>AC_stat(unterated_norm_div[i]) and treated_norm_div[i+5]>AC_stat(unterated_norm_div[i+5]):
            GCSs.append(i+1)
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
    plot_enrichment_signal(Tet_ID, treated_norm_div, unterated_norm_div, deletions, path_out)    

    #Writes GCSs data.
    GCSs_out=open(path_out+Tet_ID+'_raw_GCSs_called.txt', 'w')
    GCSs_out.write('GCSs_coordinate\tN3E\n')
    for GCS in GCSs:
        GCSs_out.write(str(GCS) + '\t' + str(min(treated_norm_div[i], treated_norm_div[i+5])) + '\n')  
    GCSs_out.close()
    return

GCSs_caller(Tetrade, Deletions, Genome, Path_for_output)

print('Script ends seccesfully!')
