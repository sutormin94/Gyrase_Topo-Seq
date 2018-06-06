###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes results of scanning procedure for forward and reverse strands.
#Returns score for every genome position (writes into WIG file), for GCSs (writes TAB files contain coordinate\tN3E\tScore info), 
#computes Pearson correlation between N3E and score, plots (Score, N3E) scatter plots.

###############################################

#######
#Packages to be imported.
#######

import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import pearsonr

#######
#Variables to be defined.
#######

#Path to the input files.
Fwd_strand_path=''
Rev_strand_path=''
#Path to the output WIG file.
Output_score_wig=''
#Dataset name for WIG header.
Dataset_name=''
#Name of the chromosome for WIG header.
Chromosome_name='NC_007779.1_w3110_Mu'
#Input data - GCSs, TAB.
path_to_GCSs_files={'Cfx': '',
                    'RifCfx': '',
                    'Micro': '',
                    'Oxo': ''}
#Output data - GCSs, TAB.
Outpath_to_GCSs_files={'Cfx': '',
                    'RifCfx': '',
                    'Micro': '',
                    'Oxo': ''}
#Output path for N3E-score plots.
N3E_score_plot_path=''

#######
#Parsing and preparing score files.
#######

def Parsing_score(path_s, path_rc, path_out, ds_name, chr_name):
    #Parse input files.
    filein_s=open(path_s, "r")
    filein_rc=open(path_rc, "r")
    fileout=open(path_out, "w")
    fileout.write("track type=wiggle_0 name=\""+str(ds_name)+"\" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom=" + str(chr_name) + " start=1 step=1\n")
    ar_s=[]
    ar_rc=[]
    for line in filein_s:
        line=line.rstrip().split('\t')
        ar_s.append([int(line[0]), float(line[1]), str(line[2])])
    for line in filein_rc:
        line=line.rstrip().split('\t')
        ar_rc.append([int(line[0]), float(line[1]), str(line[2])])    
    print(len(ar_s)) 
    print(len(ar_rc))
    #Make combined score track.
    ar_max=[]
    for i in range(4647454):
        ar_max.append(0)
    ar_max[61]=ar_s[0][1]
    for i in range(len(ar_s)-1):
        ar_max[62+i]=max(ar_s[i+1][1], ar_rc[i][1])
    ar_max[len(ar_s)+61]=ar_rc[-1][1]
    #Write score containing WIG file.
    for i in range(len(ar_max)):
        fileout.write(str(ar_max[i]) + "\n")
    filein_s.close()
    filein_rc.close()
    fileout.close()
    return ar_max

#######
#Trusted GCSs data parsing.
#Calculate and return score for GCSs.
#Write GCSs files with scores info obtained.
#######

def GCSs_parsing_score_returning_info_writing(input_dict, score_ar, output_dict):
    GCSs_sets_dict={}
    for k, v in input_dict.items():
        GCSs_dict={}
        filein=open(v, 'r')
        fileout=open(output_dict[k], 'w')
        fileout.write('GCSs_coordinate\tN3E\tScore\n')
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in ['GCSs_coordinate']:
                GCSs_dict[int(line[0])]=[float(line[1]), score_ar[int(line[0])-1]]
                fileout.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str(score_ar[int(line[0])-1]) + '\n')
            else:
                continue
        filein.close()
        fileout.close()
        GCSs_sets_dict[k]=GCSs_dict
        print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(GCSs_dict)))
    return GCSs_sets_dict

#######
#Calculate N3E-score Pearson correlation.
#######

def correlation_h_s(GCSs_sets_dict):
    HS_dict={}
    for ab, GCSs_set in GCSs_sets_dict.items():
        #[[N3E data], [Score data]]
        HS_dict[ab]=[] 
        N3E=[]
        Score=[]
        for k, v in GCSs_set.items():
            N3E.append(v[0])
            Score.append(v[1])
        #Sorting.
        Score, N3E = (list(t) for t in zip(*sorted(zip(Score, N3E))))
        print('Paerson correlation (N3E, score) for ' + str(ab) + ' : ' + str(scipy.stats.pearsonr(N3E, Score)))
        HS_dict[ab].append(N3E)
        HS_dict[ab].append(Score)
    return HS_dict

#######
#Visualize (Score, N3E) dependencies.
#######

def Plot_N3E_score(HS_dict, plot_path):
    fig=plt.figure(figsize=(15,15), dpi=300)
    #Cfx
    plot0=plt.subplot2grid((2,2),(0,0))    
    fit=np.polyfit(HS_dict['Cfx'][1], HS_dict['Cfx'][0], 1)
    fit_fn=np.poly1d(fit)  
    plot0.plot(HS_dict['Cfx'][1], HS_dict['Cfx'][0], 'o', fillstyle='none', color='#ff545a', markeredgecolor='#ff545a', markersize=2, alpha=0.8)
    plot0.plot(HS_dict['Cfx'][1], fit_fn(HS_dict['Cfx'][1]), '--k')
    plot0.set_xlabel('GCSs score', size=17) 
    plot0.set_ylabel('GCSs height', size=17)
    plot0.set_title('Cfx')
    #RifCfx
    plot1=plt.subplot2grid((2,2),(0,1))    
    fit=np.polyfit(HS_dict['RifCfx'][1], HS_dict['RifCfx'][0], 1)
    fit_fn=np.poly1d(fit) 
    plot1.plot(HS_dict['RifCfx'][1], HS_dict['RifCfx'][0], 'o', fillstyle='none', color='#fa585a', markeredgecolor='#fa585a', markersize=2, alpha=0.8)
    plot1.plot(HS_dict['RifCfx'][1], fit_fn(HS_dict['RifCfx'][1]), '--k')
    plot1.set_xlabel('GCSs score', size=17)  
    plot1.set_ylabel('GCSs height', size=17)
    plot1.set_title('RifCfx')
    #Micro
    plot2=plt.subplot2grid((2,2),(1,0))    
    fit=np.polyfit(HS_dict['Micro'][1], HS_dict['Micro'][0],1)
    fit_fn=np.poly1d(fit)   
    plot2.plot(HS_dict['Micro'][1], HS_dict['Micro'][0], 'o', fillstyle='none', color='#49ba41', markeredgecolor='#49ba41', markersize=2, alpha=0.8)
    plot2.plot(HS_dict['Micro'][1], fit_fn(HS_dict['Micro'][1]), '--k')
    plot2.set_xlabel('GCSs score', size=17) 
    plot2.set_ylabel('GCSs height', size=17)
    plot2.set_title('Micro')
    #Oxo
    plot3=plt.subplot2grid((2,2),(1,1))    
    fit=np.polyfit(HS_dict['Oxo'][1], HS_dict['Oxo'][0],1)
    fit_fn=np.poly1d(fit)   
    plot3.plot(HS_dict['Oxo'][1], HS_dict['Oxo'][0], 'o', fillstyle='none', color='#5762ff', markeredgecolor='#5762ff', markersize=2, alpha=0.8)
    plot3.plot(HS_dict['Oxo'][1], fit_fn(HS_dict['Oxo'][1]), '--k')
    plot3.set_xlabel('GCSs score', size=17)
    plot3.set_ylabel('GCSs height', size=17)
    plot3.set_title('Oxo')
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, figsize=(15, 15)) 
    plt.show()
    return

#######
#Wrapps all the functions together.
#######

def Wrapper(path_s, path_rc, path_out, ds_name, chr_name, input_dict, output_dict, plot_path):
    score_track=Parsing_score(path_s, path_rc, path_out, ds_name, chr_name)
    GCSs_info_dict=GCSs_parsing_score_returning_info_writing(input_dict, score_track, output_dict)
    Sorted_for_plot_dict=correlation_h_s(GCSs_info_dict)
    Plot_N3E_score(Sorted_for_plot_dict, plot_path)
    return
    
Wrapper(Fwd_strand_path, Rev_strand_path, Output_score_wig, Dataset_name, Chromosome_name, path_to_GCSs_files, Outpath_to_GCSs_files, N3E_score_plot_path)   

print('Script ended its work succesfully!')