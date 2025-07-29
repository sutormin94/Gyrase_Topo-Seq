###############################################
##Dmitry Sutormin, 2024##
##Topo-Seq analysis##

#Script takes SAM files as input, performs QC filtering of reads 
#relying on alignment quality and a presence of the partner: 
#only reads pairs that have a score<256 are stored.
#Than the script computes coverage depth for DNA chains separately and for both.
#Additionally it calculates N5E (number of DNA fragments starts) and 
#N3E (number of DNA fragments ends) values for every genome position strand specifically.
#Coverage depth, N3E and N5E info returns as WIG files.
#Scipt can be used for multi-chromosome referense genomes (e.g., draft bacterial genomes or eukaryotic genomes).
###############################################

#######
#Packages to be imported.
#######

import os
from os import listdir
import numpy as np
from Bio import SeqIO

#######
#Variables to be defined.
#######

# Path to the working directory
pwd="/home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy"

# Path to the input raw SAM-files
sam_path=pwd + "/SAM/"

# Path to the output/input folder with SAM-files contain proper aligned reads (score<256)
edited_sam_path=pwd + "/SAM_edited/"
if not os.path.exists(edited_sam_path):
    os.makedirs(edited_sam_path)
    
# Path to the output/input TAB files
tab_path=pwd + "/TAB/"
if not os.path.exists(tab_path):
    os.makedirs(tab_path)
    
# Path to the output/input WIG files
wig_path=pwd + "/WIG/"
if not os.path.exists(wig_path):
    os.makedirs(wig_path)
    
# Path to fasta file with referense genome (to retrieve chromosomal ID and length).
ref_gen_path="/home/niagara/Storage/D_Sutormin/Hi_C_and_Seqs_2/A_Polkhovskiy/Ref_genome/GCA_000001735.2_TAIR10.1_genomic.fna"


#######
#Reads SAM file and returns proper aligned reads (FLAG score < 256).
#######
def sam_edt(in_sam_file_path, out_sam_file_path):
    sam_input=open(in_sam_file_path, 'r')
    sam_output=open(out_sam_file_path, 'w+')
    count_tot=0 #Total number of reads
    count_in=0 #Number of innormal read pairs
    for line in sam_input:
        if line[0]=="@":
            sam_output.write(line) #Header transfer		
            continue
        count_tot+=1
        line1=line.rstrip().split('\t')			
        r=bin(int(line1[1]))
        a=int(r, 2)
        b=int('100000000000', 2)	
        c=bin(a & b)
        if int(line1[1])>=256: #Cheking the alignment quality
            count_in+=1
        else:
            sam_output.write(line) #Transfer of the proper alignments
    sam_input.close()
    sam_output.close()
    print("Total number of reads: " + str(count_tot))				
    print("Number of abnormal reads alignments: " + str(count_in))
    return

#######
#Reads and checks SAM file by counting reads that form pairs.
#######
def check_sam(sam_file_path):
    sam_file=open(sam_file_path, 'r')
    count_tot=0
    count_pair=0
    count_unpair=0
    for line in sam_file:
        if line[0]=="@": #Header cheking
            continue
        count_tot+=1
        line1=line.rstrip().split('\t')
        line2=next(sam_file)
        line22=line2.rstrip().split('\t')
        if line22[0]!= line1[0]: #Pair checking
            count_unpair+=1
        else:
            count_pair+=1
    sam_file.close()
    print("Total number of reads: " + str(count_tot))
    print("Number of read pairs: " + str(count_pair))
    print("Number of not paired reads: " + str(count_unpair))
    return

#######
#Reads SAM file and create TAB file contains left-most coordinate of 
#the read alignment and observed template length (length of the DNA fragment aligned). Keep also chromosome name.
#######
def create_tab_files(input_sam_file, output_tab_file):
    sam_file=open(input_sam_file, 'r')
    outfile=open(output_tab_file, 'w+')
    for line in sam_file:
        if line[0]=="@": #Header checking	
            continue		
        line=line.rstrip().split('\t')	
        if -1500<int(line[8])<1500: #Distance between reads within the pair has to be less than 1500 bp
            outfile.write(line[2] + "\t" + line[3] + "\t" + line[8] + "\n")
    sam_file.close()
    outfile.close()
    return

#######
#Read reference genome, get chromosomal id and length.
#######
def read_ref_genome(ref_gen_path):
    
    chr_id_len_dict={}
    
    for record in SeqIO.parse(ref_gen_path, "fasta"):
        chr_id=record.id
        chr_len=len(str(record.seq))
        
        chr_id_len_dict[chr_id]=chr_len
        print(f'Identified chromosome {chr_id}, length {chr_len} nt')
    
    return chr_id_len_dict

#######
#Reads TAB files.
#######
def Tab_pars(filein):
    
    peaks_dict={}
    num_reads=0
    
    for line in filein:
        line=line.rstrip().split('\t')
        peak=[]
        chr_id=line[0]
        peak.append(int(line[1]))
        peak.append(int(line[2]))
        
        if chr_id not in peaks_dict:
            peaks_dict[chr_id]=[peak]
            num_reads+=1
        else:
            peaks_dict[chr_id].append(peak)
            num_reads+=1
                       
    return peaks_dict, num_reads

#######
#Looks through the array, contains elements such as in TAB file (left-most coordinate + DNA fragment length) and
#constructs a new list contains read pairs.
#######
def QC_reads(ars_dict):
    
    qual_pairs_dict={}
    total_qual_reads=0
    
    for chr_id, ar in ars_dict.items():
        
        if chr_id!="*": #Filter out unmapped reads.
            
            qual_pairs=[]
            
            for i in range(len(ar)-1):
                pair=[]
                if ar[i][0]!=0 and ar[i][1]!=0 and int(ar[i][1])==-int(ar[i+1][1]): #Check that two following reads form a pair
                    pair.append(ar[i])
                    pair.append(ar[i+1])
                    qual_pairs.append(pair)
                    
            qual_pairs_dict[chr_id]=qual_pairs
            total_qual_reads+=2*len(qual_pairs)
            print(f'Number of properly paired reads for chromosome {chr_id}:  {2*len(qual_pairs)}')
    
    print(f'Total number of properly paired reads:  {total_qual_reads}')
    
    return qual_pairs_dict

#######
#Looks through the array, contains read pairs and classifies them according 
#to the orientation of the pair (forward or reverse). 
#Returns dictionary contains two lists - one for FOR pairs and one for REV pairs.
#######
def Read_strand_classif(qual_read_pairs_dict):
    
    forw={}
    rev={}
    
    total_forw=0
    total_rev=0
    
    for chr_id, ar in qual_read_pairs_dict.items():
        
        for j in range(len(ar)):
            if int(ar[j][0][1])<0: #The pair is aligned in the reverse orientation (aligned to the reverse strand)
                if chr_id not in rev:
                    rev[chr_id]=[ar[j][1]]
                else:
                    rev[chr_id].append(ar[j][1])
                    
            elif int(ar[j][0][1])>0: #The pair is aligned in the forward orientation (aligned to the forward strand)
                if chr_id not in forw:
                    forw[chr_id]=[ar[j][0]] 
                else:
                    forw[chr_id].append(ar[j][0])
                    
        total_forw+=2*len(forw[chr_id])
        total_rev+=2*len(rev[chr_id])                 
        print(f'Number of reads which DNA fragments were aligned to the forward strand for chromosome {chr_id}: {2*len(forw[chr_id])}')
        print(f'Number of reads which DNA fragments were aligned to the reverse strand for chromosome {chr_id}: {2*len(rev[chr_id])}')
        
    print(f'Total number of reads which DNA fragments were aligned to the forward strand: {total_forw}')
    print(f'Total number of reads which DNA fragments were aligned to the reverse strand: {total_rev}')
        
    return {'Forward_r' : forw, 'Reverse_r' : rev}

#######
#Looks through the array, contains read pairs classified according to 
#the orientation (aligned to the for. or to the rev. strands). 
#Returns list contains left-most and right-most coordinates of the DNA fragment aligned.
#######
def Coords(ars_dict, strand):
    
    reads_num=0
    
    ars_dict_out={}
    
    if strand=="+":
        print("Alignment to the forward strand")
        
        for chr_id, ar in ars_dict.items():
            
            ar_out=[]
            
            for k in range(len(ar)):
                pair_c=[]
                pair_c.append(int(ar[k][0]))
                pair_c.append(int(ar[k][0])+int(ar[k][1])-1)
                ar_out.append(pair_c) #Left-most and right-most coordinates of the DNA fragment aligned
            
            reads_num+=len(ar_out)  
            ars_dict_out[chr_id]=ar_out
            
    elif strand=="-":
        print("Alignment to the reverse strand")
        
        for chr_id, ar in ars_dict.items():
            
            ar_out=[]        
        
            for k in range(len(ar)):
                pair_c=[]
                pair_c.append(int(ar[k][0])-1)
                pair_c.append(int(ar[k][0])+int(ar[k][1]))
                ar_out.append(pair_c) #Left-most and right-most coordinates of the DNA fragment aligned  
            
            reads_num+=len(ar_out)    
            ars_dict_out[chr_id]=ar_out
            
    return ars_dict_out, reads_num

#######
#Calculates coverage depth for every genome position using the coordinates of 
#the DNA fragments alignment (left-most and right-most). Returns list 
#contains cov depth for every genome position.
#######
def depth_counter(coords_ars_dict, chr_id_len_dict):
    
    genome_dict={}
    
    for chr_id, chr_len in chr_id_len_dict.items():
        genome_dict[chr_id]=[0]*(chr_len+999)
        
    for chr_id, chr_len in chr_id_len_dict.items():  
        
        if chr_id in coords_ars_dict:
            coords_ar=coords_ars_dict[chr_id]
            
            for i in range(len(coords_ar)-1):

                if i%500000==0:
                    print(f'{i} read pairs processed for chromosome {chr_id}')
                    
                for k in range(coords_ar[i][1]-coords_ar[i][0]):
                    genome_dict[chr_id][coords_ar[i][0]+k]=genome_dict[chr_id][coords_ar[i][0]+k]+1             
               
        else:
            
            print(f'Chromosome with id {chr_id} was not found in the alignment data!')
        
    return genome_dict

#######
#Looks through two equal-length lists contains num values.
#Creates new list with pairwise summs.
#######
def Integrator(ars_dict1, ars_dict2, chr_id_len_dict):
    
    int_genome_dict={}
    
    for chr_id, chr_len in chr_id_len_dict.items():
        int_genome_dict[chr_id]=[0]*(chr_len+999)
        
    for chr_id, chr_len in chr_id_len_dict.items():  
        
        if (chr_id in ars_dict1) and (chr_id in ars_dict2):
            ar1=ars_dict1[chr_id]
            ar2=ars_dict2[chr_id]
            
            for i in range (len(ar1)-1):
                int_genome_dict[chr_id][i]=ar1[i]+ar2[i]            
            
        else:
            
            print(f'Chromosome with id {chr_id} was not found in the alignment data!')

    return int_genome_dict

#######
#Writes WIG file using the array of ints or floats.
#######
def write_file(genome_dict, flag, chr_id_len_dict, fileout_path):
    fileout=open(fileout_path, 'w+')
    
    for chr_id, ar in genome_dict.items():
    
        fileout.write('track type=wiggle_0 name="'+str(flag)+'" autoScale=off viewLimits=0.0:25.0'+'\n'+'fixedStep chrom='+str(chr_id)+' start=1 step=1\n')   
        for i in range(len(ar)):
            fileout.write(str(ar[i])+'\n')
            
    fileout.close()
    
    return

#######
#Calculates the number of DNA fragments starts (N5E) and ends (N3E) depends on the fragment orientation 
#for every genome position. And integrates the values if start and end can not be distinguish from each other
#Returns this information as dict of lists.
#######
def start_end_count(forw_dict, rev_dict, chr_id_len_dict): 
    
    genome_start_F_dict={}
    genome_end_F_dict={}
    genome_start_R_dict={}
    genome_end_R_dict={}
    
    for chr_id, chr_len in chr_id_len_dict.items():

        genome_start_F=[0]*(chr_len+999)
        genome_end_F=[0]*(chr_len+999)
        genome_start_R=[0]*(chr_len+999)
        genome_end_R=[0]*(chr_len+999)
        
        if (chr_id in forw_dict) and (chr_id in rev_dict):
            
            forw=forw_dict[chr_id]
            rev=rev_dict[chr_id]
            
            for i in range(len(forw)):
                genome_start_F[forw[i][0]-1]=genome_start_F[forw[i][0]-1]+1
                genome_end_F[forw[i][1]-1]=genome_end_F[forw[i][1]-1]+1
            for i in range(len(rev)):
                genome_start_R[rev[i][1]]=genome_start_R[rev[i][1]]+1
                genome_end_R[rev[i][0]]=genome_end_R[rev[i][0]]+1 
                
            genome_start_F_dict[chr_id]=genome_start_F   
            genome_end_F_dict[chr_id]=genome_end_F
            genome_start_R_dict[chr_id]=genome_start_R
            genome_end_R_dict[chr_id]=genome_end_R
        
        else:
            
            print(f'Chromosome with id {chr_id} was not found in the alignment data!')
                
    genome_start_and_end_F_dict=Integrator(genome_start_F_dict, genome_end_F_dict, chr_id_len_dict)
    genome_start_and_end_R_dict=Integrator(genome_start_R_dict, genome_end_R_dict, chr_id_len_dict)
    
    return_dict={'DNA_fragments_starts_forward_strand' : genome_start_F_dict, 
                     'DNA_fragments_ends_forward_strand' : genome_end_F_dict, 
                     'DNA_fragments_starts_and_ends_forward_strand' : genome_start_and_end_F_dict,
                     'DNA_fragments_starts_reverse_strand' : genome_start_R_dict, 
                     'DNA_fragments_ends_reverse_strand' : genome_end_R_dict, 
                     'DNA_fragments_starts_and_ends_reverse_strand' : genome_start_and_end_R_dict}
            
    return return_dict

#######
#Wraps functions that read, edit and write SAM files (sam_edt) and check the resulting edited SAM (check_sam). 
#Editing results in filtering of the proper aligned reads with score<256.
#Checking procedure: counting reads that form pairs.
#######
def edit_sam_files_wrapper(sam_path, edited_sam_path, check_option):
    #Prepares the list of .sam files to work with
    files=listdir(sam_path)
    input_samfiles=[]
    for file in files:
        if file.endswith(".sam"):
            input_samfiles.append(file)	
    #Edit .sam files
    for sam in input_samfiles:
        print(sam)
        in_sam_file=sam_path + sam
        out_sam_file=edited_sam_path + sam[:-4] + "_edt.sam"
        sam_edt(in_sam_file, out_sam_file)
    #Check .sam files (optional)
    if check_option==1:
        files=listdir(edited_sam_path)
        edited_samfiles=[]
        for file in files:
            if file.endswith("_edt.sam"):
                edited_samfiles.append(file)
        for sam in edited_samfiles:
            sam_file_path=edited_sam_path + sam
            check_sam(sam_file_path)
    return

#######
#Wraps functions that read SAM files and makes TAB files (create_tab_files).
#While running, it filters reads pairs consist of reads that form a DNA fragment
#less than 1500 bp.
#######	
def create_tab_files_wrapper(edited_sam_path, tab_path):
    #Reads SAM files were edited
    files=listdir(edited_sam_path)
    samfiles=[]
    for file in files:
        if file.endswith(".sam"):
            samfiles.append(file)		
    #Creates TAB files
    for sam in samfiles:
        print(sam)
        input_sam_file_path=edited_sam_path + sam
        out_sam_file_path=tab_path + sam[:-4] + ".tab"
        create_tab_files(input_sam_file_path, out_sam_file_path)
    return

#######
#Wraps functions that read TAB files (Peaks_pars from pars_com),
#makes reads pairs (QC_reads), strand classify reads pairs (Read_strand_classif),
#marks left- and right-most positions of the DNA fragments aligned (Coords),
#calculates coverage depth for + and - strands and sum them (depth_counter, Integrator), 
#calculates number of DNA fragments starts (N5E) and ends (N3E) for every genome position (start_end_count) and
#writes output WIG files
#######
def create_wig_files_wrapper(tab_path, wig_path, ref_gen_path):
    
    #Read referense genome, get chromosomal id and length.
    chr_id_len_dict=read_ref_genome(ref_gen_path)
    #Makes the list of TAB files to work with
    files=listdir(tab_path)
    tabfiles=[]
    for file in files:
        if file.endswith(".tab"):
            tabfiles.append(file)
    #TAB files parsing, pairs construction, strand classification, start-end calculation,
    #coverage depth calculation, coverage depth integration, number of start-end calculation,
    #WIG files writing.
    for tab in tabfiles:
        #Parsing of the .tab file
        print(tab)
        tab_file=open(tab_path+tab, 'r')
        reads_ars_dict, num_reads=Tab_pars(tab_file)
        print("Number of reads in the " + str(tab) + " file: " + str(num_reads))
        #Reads pairs
        qual_read_pairs_dict=QC_reads(reads_ars_dict)
        #Read pairs classification according to the alignment orientation
        classified_reads=Read_strand_classif(qual_read_pairs_dict) 
        Forward_dict=classified_reads['Forward_r']
        Reverse_dict=classified_reads['Reverse_r']
        #Left-most and right-most positions of the DNA fragments aligned
        for_coords_dict, reads_num_for=Coords(Forward_dict, "+")
        rev_coords_dict, reads_num_rev=Coords(Reverse_dict, "-")
        print("Number of paired reads for file " + str(tab) + "\nthat forms DNA fragments aligned to the forward strand: " + str(2*reads_num_for))
        print("Number of paired reads for file " + str(tab) + "\nthat forms DNA fragments aligned to the reverse strand: " + str(2*reads_num_rev))
        #Calculates coverage depth for 2 strands separately
        for_reads_genome_depth_dict=depth_counter(for_coords_dict, chr_id_len_dict)  
        rev_reads_genome_depth_dict=depth_counter(rev_coords_dict, chr_id_len_dict)
        #Calculates sum coverage depth for both strands
        genome_depth_dict=Integrator(for_reads_genome_depth_dict, rev_reads_genome_depth_dict, chr_id_len_dict)

        #Writes WIG files (coverage depth)
        outfile_path=wig_path + tab[:-4] + "_forward_depth.wig"
        write_file(for_reads_genome_depth_dict, str(tab[:-4]) + "for_reads", chr_id_len_dict, outfile_path)
        outfile_path=wig_path + tab[:-4] + "_reverse_depth.wig"
        write_file(rev_reads_genome_depth_dict, str(tab[:-4]) + "rev_reads", chr_id_len_dict, outfile_path)
        outfile_path=wig_path + tab[:-4] + "_for_rev_depth.wig"
        write_file(genome_depth_dict, str(tab[:-4]) + "depth", chr_id_len_dict, outfile_path)

        #Calculates number of DNA fragments starts (N5E) and ends (N3E) for every genome position
        N3E_N5E_dict=start_end_count(for_coords_dict, rev_coords_dict, chr_id_len_dict)		

        Starts_F_dict=N3E_N5E_dict['DNA_fragments_starts_forward_strand']
        Ends_F_dict=N3E_N5E_dict['DNA_fragments_ends_forward_strand']
        Starts_and_ends_F_dict=N3E_N5E_dict['DNA_fragments_starts_and_ends_forward_strand']
        Starts_R_dict=N3E_N5E_dict['DNA_fragments_starts_reverse_strand']
        Ends_R_dict=N3E_N5E_dict['DNA_fragments_ends_reverse_strand']
        Starts_and_ends_R_dict=N3E_N5E_dict['DNA_fragments_starts_and_ends_reverse_strand']		

        #Writes WIG files (N3E and N5E)
        outfile_starts_F_path=wig_path + tab[:-4] + "_N5E_F.wig" #N5E forward strand
        write_file(Starts_F_dict, str(tab[:-4]) + "_N5E_F", chr_id_len_dict, outfile_starts_F_path) 
        outfile_ends_F_path=wig_path + tab[:-4] + "_N3E_F.wig" #N3E forward strand
        write_file(Ends_F_dict, str(tab[:-4]) + "_N3E_F", chr_id_len_dict, outfile_ends_F_path) 
        outfile_starts_ends_F_path=wig_path + tab[:-4] + "_N53E_F.wig" #N5E+N3E forward strand
        write_file(Starts_and_ends_F_dict, str(tab[:-4]) + "_N5E_and_N3E_F", chr_id_len_dict, outfile_starts_ends_F_path) 

        outfile_starts_R_path=wig_path + tab[:-4] + "_N5E_R.wig" #N5E reverse strand
        write_file(Starts_R_dict, str(tab[:-4]) + "_N5E_R", chr_id_len_dict, outfile_starts_R_path) 
        outfile_ends_R_path=wig_path + tab[:-4] + "_N3E_R.wig" #N3E reverse strand
        write_file(Ends_R_dict, str(tab[:-4]) + "_N3E_R", chr_id_len_dict, outfile_ends_R_path) 
        outfile_starts_ends_R_path=wig_path + tab[:-4] + "_N53E_R.wig" #N5E+N3E reverse strand
        write_file(Starts_and_ends_R_dict, str(tab[:-4]) + "_N5E_and_N3E_R", chr_id_len_dict, outfile_starts_ends_R_path) 		

    return


edit_sam_files_wrapper(sam_path, edited_sam_path, 0)
create_tab_files_wrapper(edited_sam_path, tab_path)
create_wig_files_wrapper(tab_path, wig_path, ref_gen_path)

print('Script ended its work succesfully!')
