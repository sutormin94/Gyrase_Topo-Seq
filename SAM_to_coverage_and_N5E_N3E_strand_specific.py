###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Script takes SAM files as input, performs QC filtering of reads 
#relying on alignment quality and a presence of the partner: 
#only reads pairs that have a score<256 are stored.
#Than the script computes coverage depth for DNA chains separately and for both.
#Additionally it calculates N5E (number of DNA fragments starts) and 
#N3E (number of DNA fragments ends) values for every genome position.
#Coverage depth, N3E and N5E info returns as WIG files.
###############################################

#######
#Packages to be imported.
#######

import os
from os import listdir
import numpy as np

#######
#Variables to be defined.
#######

#Path to the working directory
pwd="/data/Gyrase/Data_preparation/Cfx_900nM"
#Path to the input raw SAM-files
sam_path=pwd + "/SAM/"
#Path to the output/input folder with SAM-files contain proper aligned reads (score<256)
edited_sam_path=pwd + "/SAM_edited/"
if not os.path.exists(edited_sam_path):
	os.makedirs(edited_sam_path)
#Path to the output/input TAB files
tab_path=pwd + "/TAB/"
if not os.path.exists(tab_path):
	os.makedirs(tab_path)
#Path to the output/input WIG files
wig_path=pwd + "/WIG/"
if not os.path.exists(wig_path):
	os.makedirs(wig_path)
#Chromosome (genome) identificator (look for the corresponding FASTA ID)
chromosome_identificator="NC_007779.1_w3110_Mu"
chromosome_length=4647999


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
		c=bin (a & b)
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
#the read alignment and observed template length (length of the DNA fragment aligned).
#######
def create_tab_files(input_sam_file, output_tab_file):
	sam_file=open(input_sam_file, 'r')
	outfile=open(output_tab_file, 'w+')
	for line in sam_file:
		if line[0]=="@": #Header checking	
			continue		
		line=line.rstrip().split('\t')	
		if -1500<int(line[8])<1500: #Distance between reads within the pair has to be less than 1500 bp
			outfile.write(line[3] + "\t" + line[8] + "\n") 
	sam_file.close()
	outfile.close()
	return

#######
#Reads TAB files.
#######
def Tab_pars(filein):
	peaks=[]
	for line in filein:
		line=line.rstrip().split('\t')
		peak=[]
		peak.append(int(line[0]))
		peak.append(int(line[1]))
		peaks.append(peak)            
	return peaks

#######
#Looks through the array, contains elements such as in TAB file (left-most coordinate + DNA fragment length) and
#constructs a new list contains read pairs.
#######
def QC_reads(ar):
	qual_pairs=[]
	for i in range(len(ar)-1):
		pair=[]
		if ar[i][0]!=0 and ar[i][1]!=0 and int(ar[i][1])==-int(ar[i+1][1]): #Check that two following reads form a pair
			pair.append(ar[i])
			pair.append(ar[i+1])
			qual_pairs.append(pair)
	print("Number of proper paired reads: " + str(2*len(qual_pairs)))
	return qual_pairs

#######
#Looks through the array, contains read pairs and classifies them according 
#to the orientation of the pair (forward or reverse). 
#Returns dictionary contains two lists - one for FOR pairs and one for REV pairs.
#######
def Read_strand_classif(ar):
	forw=[]
	rev=[]
	for j in range(len(ar)):
		if int(ar[j][0][1])<0: #The pair is aligned in the reverse orientation (aligned to the reverse strand)
			rev.append(ar[j][1])
		elif int(ar[j][0][1])>0: #The pair is aligned in the forward orientation (aligned to the forward strand)
			forw.append(ar[j][0])
	print("Number of reads which DNA fragments were aligned to the forward strand: " + str(2*len(forw)))
	print("Number of reads which DNA fragments were aligned to the reverse strand: " + str(2*len(rev)))
	return {'Forward_r':forw, 'Reverse_r':rev}

#######
#Looks through the array, contains read pairs classified according to 
#the orientation (aligned to the for. or to the rev. strands). 
#Returns list contains left-most and right-most coordinates of the DNA fragment aligned.
#######
def Coords(ar, strand):
	ar_out=[]
	if strand=="+":
		print("Alignment to the forward strand")
		for k in range(len(ar)):
			pair_c=[]
			pair_c.append(int(ar[k][0]))
			pair_c.append(int(ar[k][0])+int(ar[k][1])-1)
			ar_out.append(pair_c) #Left-most and right-most coordinates of the DNA fragment aligned
	elif strand=="-":
		print("Alignment to the reverse strand")
		for k in range(len(ar)):
			pair_c=[]
			pair_c.append(int(ar[k][0])-1)
			pair_c.append(int(ar[k][0])+int(ar[k][1]))
			ar_out.append(pair_c) #Left-most and right-most coordinates of the DNA fragment aligned    
	return ar_out

#######
#Calculates coverage depth for every genome position using the coordinates of 
#the DNA fragments alignment (left-most and right-most). Returns list 
#contains cov depth for every genome position.
#######
def depth_counter(coords_ar, chromosome_length):
	genome=[]
	for i in range(chromosome_length):
		genome.append(0)
	i=0
	counter=np.arange(0, 8000000, 100000)
	for i in range(len(coords_ar)-1):
		i=i+1
		if i in counter:
			print(i)
		for k in range(coords_ar[i][1]-coords_ar[i][0]):
			genome[coords_ar[i][0]+k]=genome[coords_ar[i][0]+k]+1     
	return genome

#######
#Looks through two equal-length lists contains num values.
#Creates new list with pairwise summs.
#######
def Integrator(ar1, ar2, chromosome_length):
	ar3=[0]*chromosome_length
  
	for i in range (len(ar1)-1):
		ar3[i]=ar1[i]+ar2[i]
	return ar3

#######
#Writes WIG file using the array of ints or floats.
#######
def write_file(ar, flag, strain_id, fileout_path):
	fileout=open(fileout_path, 'w+')
	fileout.write('track type=wiggle_0 name="'+str(flag)+'" autoScale=off viewLimits=0.0:25.0'+'\n'+'fixedStep chrom='+str(strain_id)+' start=1 step=1\n')   
	for i in range(len(ar)):
		fileout.write(str(ar[i])+'\n')
	fileout.close()
	return

#######
#Calculates the number of DNA fragments starts (N5E) and ends (N3E) depends on the fragment orientation 
#for every genome position. And integrates the values if start and end can not be distinguish from each other
#Returns this information as dict of lists.
#######
def start_end_count(forw, rev, chromosome_length):  
	
	genome_start_F=[0]*chromosome_length
	genome_end_F=[0]*chromosome_length
	genome_start_R=[0]*chromosome_length
	genome_end_R=[0]*chromosome_length
	
	for i in range(len(forw)):
		genome_start_F[forw[i][0]-1]=genome_start_F[forw[i][0]-1]+1
		genome_end_F[forw[i][1]-1]=genome_end_F[forw[i][1]-1]+1
	for i in range(len(rev)):
		genome_start_R[rev[i][1]]=genome_start_R[rev[i][1]]+1
		genome_end_R[rev[i][0]]=genome_end_R[rev[i][0]]+1 
	genome_start_and_end_F=Integrator(genome_start_F, genome_end_F, chromosome_length)
	genome_start_and_end_R=Integrator(genome_start_R, genome_end_R, chromosome_length)
	
	return_dict={'DNA_fragments_starts_forward_strand' : genome_start_F, 
		     'DNA_fragments_ends_forward_strand' : genome_end_F, 
		     'DNA_fragments_starts_and_ends_forward_strand' : genome_start_and_end_F,
		     'DNA_fragments_starts_reverse_strand' : genome_start_R, 
		     'DNA_fragments_ends_reverse_strand' : genome_end_R, 
		     'DNA_fragments_starts_and_ends_reverse_strand' : genome_start_and_end_R}
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
def create_wig_files_wrapper(tab_path, wig_path, chromosome_id, chromosome_length):
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
		reads_ar=Tab_pars(tab_file)
		print("Number of reads in the " + str(tab) + " file: " + str(len(reads_ar)))
		#Reads pairs
		qual_read_pairs=QC_reads(reads_ar)
		#Read pairs classification according to the alignment orientation
		classified_reads=Read_strand_classif(qual_read_pairs) 
		Forward=classified_reads['Forward_r']
		Reverse=classified_reads['Reverse_r']
		#Left-most and right-most positions of the DNA fragments aligned
		for_coords=Coords(Forward, "+")
		rev_coords=Coords(Reverse, "-")
		print("Number of paired reads for file " + str(tab) + "\nthat forms DNA fragments aligned to the forward strand: " + str(2*len(for_coords)))
		print("Number of paired reads for file " + str(tab) + "\nthat forms DNA fragments aligned to the reverse strand: " + str(2*len(rev_coords)))
		#Calculates coverage depth for 2 strands separately
		for_reads_genome_depth=depth_counter(for_coords, chromosome_length)  
		rev_reads_genome_depth=depth_counter(rev_coords, chromosome_length)
		#Calculates sum coverage depth for both strands
		genome_depth=Integrator(for_reads_genome_depth, rev_reads_genome_depth, chromosome_length)
		
		#Writes WIG files (coverage depth)
		outfile_path=wig_path + tab[:-4] + "_forward_depth.wig"
		write_file(for_reads_genome_depth, str(tab[:-4]) + "for_reads", chromosome_id, outfile_path)
		outfile_path=wig_path + tab[:-4] + "_reverse_depth.wig"
		write_file(rev_reads_genome_depth, str(tab[:-4]) + "rev_reads", chromosome_id, outfile_path)
		outfile_path=wig_path + tab[:-4] + "_for_rev_depth.wig"
		write_file(genome_depth, str(tab[:-4]) + "depth", chromosome_id, outfile_path)
		
		#Calculates number of DNA fragments starts (N5E) and ends (N3E) for every genome position
		N3E_N5E_dict=start_end_count(for_coords, rev_coords, chromosome_length)		
		
		Starts_F=N3E_N5E_dict['DNA_fragments_starts_forward_strand']
		Ends_F=N3E_N5E_dict['DNA_fragments_ends_forward_strand']
		Starts_and_ends_F=N3E_N5E_dict['DNA_fragments_starts_and_ends_forward_strand']
		Starts_R=N3E_N5E_dict['DNA_fragments_starts_reverse_strand']
		Ends_R=N3E_N5E_dict['DNA_fragments_ends_reverse_strand']
		Starts_and_ends_R=N3E_N5E_dict['DNA_fragments_starts_and_ends_reverse_strand']		
		
		#Writes WIG files (N3E and N5E)
		outfile_starts_F_path=wig_path + tab[:-4] + "_N5E_F.wig" #N5E forward strand
		write_file(Starts_F, str(tab[:-4]) + "_N5E_F", chromosome_id, outfile_starts_F_path) 
		outfile_ends_F_path=wig_path + tab[:-4] + "_N3E_F.wig" #N3E forward strand
		write_file(Ends_F, str(tab[:-4]) + "_N3E_F", chromosome_id, outfile_ends_F_path) 
		outfile_starts_ends_F_path=wig_path + tab[:-4] + "_N53E_F.wig" #N5E+N3E forward strand
		write_file(Starts_and_ends_F, str(tab[:-4]) + "_N5E_and_N3E_F", chromosome_id, outfile_starts_ends_F_path) 
		
		outfile_starts_R_path=wig_path + tab[:-4] + "_N5E_R.wig" #N5E reverse strand
		write_file(Starts_R, str(tab[:-4]) + "_N5E_R", chromosome_id, outfile_starts_R_path) 
		outfile_ends_R_path=wig_path + tab[:-4] + "_N3E_R.wig" #N3E reverse strand
		write_file(Ends_R, str(tab[:-4]) + "_N3E_R", chromosome_id, outfile_ends_R_path) 
		outfile_starts_ends_R_path=wig_path + tab[:-4] + "_N53E_R.wig" #N5E+N3E reverse strand
		write_file(Starts_and_ends_R, str(tab[:-4]) + "_N5E_and_N3E_R", chromosome_id, outfile_starts_ends_R_path) 		
			
	return


edit_sam_files_wrapper(sam_path, edited_sam_path, 0)
create_tab_files_wrapper(edited_sam_path, tab_path)
create_wig_files_wrapper(tab_path, wig_path, chromosome_identificator, chromosome_length)

print('Script ended its work succesfully!')
