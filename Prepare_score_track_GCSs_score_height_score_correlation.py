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

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import pearsonr

#######
#Variables to be defined.
#######

#Path to the input files.
Fwd_strand_path="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\Combined_motif\E_coli_w3110_G_Mu_scan_forward_with_combined_motif.txt"
Rev_strand_path="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\Combined_motif\E_coli_w3110_G_Mu_scan_rc_with_combined_motif.txt"

#Path to the output WIG file.
path_to_res_score_files="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\Score_tracks\\"
if not os.path.exists(path_to_res_score_files):
    os.makedirs(path_to_res_score_files)
Output_score_wig=path_to_res_score_files + "E_coli_w3110_G_Mu_score.wig"
#Dataset name for WIG header.
Dataset_name="Cfx_Micro_Oxo_1828"
#Name of the chromosome for WIG header.
Chromosome_name='NC_007779.1_w3110_Mu'

#Input data - GCSs, TAB.
path_to_GCSs_files={'Cfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\Cfx_10mkM_trusted_GCSs.txt",
                    'RifCfx': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\RifCfx_trusted_GCSs.txt",
                    'Micro': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\Micro_trusted_GCSs.txt",
                    'Oxo': "C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets\Oxo_trusted_GCSs.txt"}
#Output data - GCSs, TAB.
path_to_GCSs_sets_with_score="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\GCSs_sets_score\\"
if not os.path.exists(path_to_GCSs_sets_with_score):
    os.makedirs(path_to_GCSs_sets_with_score)
Outpath_to_GCSs_files={'Cfx': path_to_GCSs_sets_with_score + "Cfx_10mkM_trusted_GCSs_h_s.txt",
                    'RifCfx': path_to_GCSs_sets_with_score + "RifCfx_trusted_GCSs_h_s.txt",
                    'Micro': path_to_GCSs_sets_with_score + "Micro_trusted_GCSs_h_s.txt",
                    'Oxo': path_to_GCSs_sets_with_score + "Oxo_trusted_GCSs_h_s.txt"}
#Output directory for N3E-score plots.
N3E_score_plot_path=path_to_GCSs_sets_with_score

#######
#Parsing and preparing score files.
#######

def Parsing_score(path_s, path_rc, path_out, ds_name, chr_name):
    #Parse input files.
    filein_s=open(path_s, "r")
    filein_rc=open(path_rc, "r")
    fileout=open(path_out, "w")
    fileout.write("track type=wiggle_0 name=\""+str(ds_name)+"\" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom=" + str(chr_name) + " start=1 step=1\n")
    ar_s={}
    ar_rc={}
    for line in filein_s:
        line=line.rstrip().split('\t')
        ar_s[int(line[0])]=[float(line[1]), str(line[2])]
    for line in filein_rc:
        line=line.rstrip().split('\t')
        ar_rc[int(line[0])]=[float(line[1]), str(line[2])]    
    print("Length of the fw scan track: " + str(len(ar_s)))
    print("Length of the rc scan track: " + str(len(ar_rc)))
    #Make combined score track.
    ar_max=[]
    for i in range(4647454):
        if i+1 in ar_s and i+1 in ar_rc:
            ar_max.append(max(ar_s[i+1][0], ar_rc[i+1][0]))
        elif i+1 in ar_s and i+1 not in ar_rc:
            ar_max.append(ar_s[i+1][0])
        elif i+1 not in ar_s and i+1 in ar_rc:
            ar_max.append(ar_rc[i+1][0])
        elif i+1 not in ar_s and i+1 not in ar_rc:
            ar_max.append(0)
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
                GCSs_dict[int(line[0])]=[float(line[1]), (score_ar[int(line[0])-1+1]+score_ar[int(line[0])-1+4])/2]
                fileout.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str((score_ar[int(line[0])-1+1]+score_ar[int(line[0])-1+4])/2) + '\n')
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
    fig=plt.figure(figsize=(15,15), dpi=100)
    #Cfx
    plot0=plt.subplot2grid((2,2),(0,0))    
    fit=np.polyfit(HS_dict['Cfx'][1], HS_dict['Cfx'][0], 1)
    print(fit)
    fit_fn=np.poly1d(fit)  
    plot0.plot(HS_dict['Cfx'][1], HS_dict['Cfx'][0], 'o', fillstyle='none', color='#7FCE79', markeredgecolor='#7FCE79', markersize=2, alpha=0.8)
    plot0.plot(HS_dict['Cfx'][1], fit_fn(HS_dict['Cfx'][1]), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot0.set_xlabel('GCSs score', size=17) 
    plot0.set_ylabel('GCSs N3E', size=17)
    plot0.legend(loc='upper right')
    plot0.set_title('Cfx', size=18)
    #RifCfx
    plot1=plt.subplot2grid((2,2),(0,1))    
    fit=np.polyfit(HS_dict['RifCfx'][1], HS_dict['RifCfx'][0], 1)
    print(fit)
    fit_fn=np.poly1d(fit) 
    plot1.plot(HS_dict['RifCfx'][1], HS_dict['RifCfx'][0], 'o', fillstyle='none', color='#BAE85C', markeredgecolor='#BAE85C', markersize=2, alpha=0.8)
    plot1.plot(HS_dict['RifCfx'][1], fit_fn(HS_dict['RifCfx'][1]), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot1.set_xlabel('GCSs score', size=17)  
    plot1.set_ylabel('GCSs N3E', size=17)
    plot1.legend(loc='upper right')
    plot1.set_title('RifCfx', size=18)
    #Micro
    plot2=plt.subplot2grid((2,2),(1,0))    
    fit=np.polyfit(HS_dict['Micro'][1], HS_dict['Micro'][0],1)
    print(fit)
    fit_fn=np.poly1d(fit)   
    plot2.plot(HS_dict['Micro'][1], HS_dict['Micro'][0], 'o', fillstyle='none', color='#ff878b', markeredgecolor='#ff878b', markersize=2, alpha=0.8)
    plot2.plot(HS_dict['Micro'][1], fit_fn(HS_dict['Micro'][1]), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot2.set_xlabel('GCSs score', size=17) 
    plot2.set_ylabel('GCSs N3E', size=17)
    plot2.legend(loc='upper right')
    plot2.set_title('Micro', size=18)
    #Oxo
    plot3=plt.subplot2grid((2,2),(1,1))    
    fit=np.polyfit(HS_dict['Oxo'][1], HS_dict['Oxo'][0],1)
    print(fit)
    fit_fn=np.poly1d(fit)   
    plot3.plot(HS_dict['Oxo'][1], HS_dict['Oxo'][0], 'o', fillstyle='none', color='#8991ff', markeredgecolor='#8991ff', markersize=2, alpha=0.8)
    plot3.plot(HS_dict['Oxo'][1], fit_fn(HS_dict['Oxo'][1]), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
    plot3.set_xlabel('GCSs score', size=17)
    plot3.set_ylabel('GCSs N3E', size=17)
    plot3.legend(loc='upper right')
    plot3.set_title('Oxo', size=18)
    plt.tight_layout()
    plt.savefig(plot_path + "GCSs_N3E_score_correlations.png", dpi=400, figsize=(15, 15)) 
    plt.close()
    return

#######
#Visualize Score, N3E distributions.
#######

def Plot_N3E_score_distribution(HS_dict, score_ar, plot_path):
    #matplotlib.rc('text', usetex = True)
    fig=plt.figure(figsize=(15,15), dpi=100)
    #Cfx N3E
    plot0=plt.subplot2grid((4,3),(0,0), rowspan=1, colspan=1)     
    plot0.hist(HS_dict['Cfx'][0], color='#7FCE79', edgecolor='black', alpha=0.8, label='Mean N3E='+str(round(np.mean(HS_dict['Cfx'][0]),2)))
    plot0.set_xlabel('GCSs N3E', size=17) 
    plot0.set_ylabel('Number of GCSs', size=17)
    plot0.legend(loc='upper right')
    plot0.set_title('Cfx', size=18)
    #Cfx score
    plot1=plt.subplot2grid((4,3),(1,0), rowspan=1, colspan=1)     
    plot1.hist(HS_dict['Cfx'][1], color='#7FCE79', edgecolor='black', alpha=0.8, label='Mean score='+str(round(np.mean(HS_dict['Cfx'][1]),2)))
    plot1.set_xlabel('GCSs score', size=17) 
    plot1.set_ylabel('Number of GCSs', size=17)
    plot1.legend(loc='upper right')
    #RifCfx N3E
    plot2=plt.subplot2grid((4,3),(0,1), rowspan=1, colspan=1)     
    plot2.hist(HS_dict['RifCfx'][0], color='#BAE85C', edgecolor='black', alpha=0.8, label='Mean N3E='+str(round(np.mean(HS_dict['RifCfx'][0]),2)))
    plot2.set_xlabel('GCSs N3E', size=17) 
    plot2.set_ylabel('Number of GCSs', size=17)
    plot2.legend(loc='upper right')
    plot2.set_title('RifCfx', size=18)
    #RifCfx score
    plot3=plt.subplot2grid((4,3),(1,1), rowspan=1, colspan=1)     
    plot3.hist(HS_dict['RifCfx'][1], color='#BAE85C', edgecolor='black', alpha=0.8, label='Mean score='+str(round(np.mean(HS_dict['RifCfx'][1]),2)))
    plot3.set_xlabel('GCSs score', size=17) 
    plot3.set_ylabel('Number of GCSs', size=17)  
    plot3.legend(loc='upper right')
    #Micro N3E
    plot4=plt.subplot2grid((4,3),(0,2), rowspan=1, colspan=1)     
    plot4.hist(HS_dict['Micro'][0], color='#ff878b', edgecolor='black', alpha=0.8, label='Mean N3E='+str(round(np.mean(HS_dict['Micro'][0]),2)))
    plot4.set_xlabel('GCSs N3E', size=17) 
    plot4.set_ylabel('Number of GCSs', size=17)
    plot4.legend(loc='upper right')
    plot4.set_title('Micro', size=18)
    #Micro score
    plot5=plt.subplot2grid((4,3),(1,2), rowspan=1, colspan=1)     
    plot5.hist(HS_dict['Micro'][1], color='#ff878b', edgecolor='black', alpha=0.8, label='Mean score='+str(round(np.mean(HS_dict['Micro'][1]),2)))
    plot5.set_xlabel('GCSs score', size=17) 
    plot5.set_ylabel('Number of GCSs', size=17)  
    plot5.legend(loc='upper right')
    #Oxo N3E
    plot6=plt.subplot2grid((4,3),(2,0), rowspan=1, colspan=1)     
    plot6.hist(HS_dict['Oxo'][0], color='#8991ff', edgecolor='black', alpha=0.8, label='Mean N3E='+str(round(np.mean(HS_dict['Oxo'][0]),2)))
    plot6.set_xlabel('GCSs N3E', size=17) 
    plot6.set_ylabel('Number of GCSs', size=17)
    plot6.legend(loc='upper right')
    plot6.set_title('Oxo', size=18)
    #Oxo score
    plot7=plt.subplot2grid((4,3),(3,0), rowspan=1, colspan=1)     
    plot7.hist(HS_dict['Oxo'][1], color='#8991ff', edgecolor='black', alpha=0.8, label='Mean score='+str(round(np.mean(HS_dict['Oxo'][1]),2)))
    plot7.set_xlabel('GCSs score', size=17) 
    plot7.set_ylabel('Number of GCSs', size=17)
    plot7.legend(loc='upper right')
    #Whole genome scores distribution
    plot8=plt.subplot2grid((4,3),(2,1), rowspan=2, colspan=2)  
    plot8.hist(score_ar, color='#5762ff', edgecolor='black', alpha=0.8, label='Mean score='+str(round(np.mean(score_ar),2)))
    plot8.set_xlabel('Score', size=17) 
    plot8.set_ylabel('Number of genome positions', size=17)
    plot8.legend(loc='upper right')
    plot8.set_title('Scores distribution for $\it{E. coli}$ DY330 MuSGS genome', size=18)
    plt.tight_layout()
    plt.savefig(plot_path + "GCSs_N3E_score_distributions.png", dpi=400, figsize=(15, 15)) 
    plt.close()    
    return

#######
#Wrapps all the functions together.
#######

def Wrapper(path_s, path_rc, path_out, ds_name, chr_name, input_dict, output_dict, plot_path):
    score_track=Parsing_score(path_s, path_rc, path_out, ds_name, chr_name)
    GCSs_info_dict=GCSs_parsing_score_returning_info_writing(input_dict, score_track, output_dict)
    Sorted_for_plot_dict=correlation_h_s(GCSs_info_dict)
    Plot_N3E_score(Sorted_for_plot_dict, plot_path)
    Plot_N3E_score_distribution(Sorted_for_plot_dict, score_track, plot_path)
    return
    
Wrapper(Fwd_strand_path, Rev_strand_path, Output_score_wig, Dataset_name, Chromosome_name, path_to_GCSs_files, Outpath_to_GCSs_files, N3E_score_plot_path)   

print('Script ended its work succesfully!')