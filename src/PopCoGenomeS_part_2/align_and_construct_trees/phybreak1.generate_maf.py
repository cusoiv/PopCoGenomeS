import os
import sys
import numpy

## Collect parameters
project_dir = ""
input_contig_dir = ""
contig_dir = ""
contig_extension = ""
output_prefix = ""
pop_infile_name = ""
ref_iso = ""
ref_contig = ""
len_block_threshold = 0
gap_prop_thresh = 0.0
MUGSY_source = ""


parameter_file = open("phybreak_parameters.txt","r")
for line in parameter_file:
	line = line.strip().split(" = ")
	if len(line) > 1:
		if line[0] == "project_dir":
			project_dir = line[1].split(" #")[0]
		elif line[0] == "input_contig_dir":
			input_contig_dir = line[1].split(" #")[0]
		elif line[0] == "contig_dir":
			contig_dir = line[1].split(" #")[0]
		elif line[0] == "input_contig_extension":
			contig_extension = line[1].split(" #")[0]
		elif line[0] == "output_prefix":
			output_prefix = line[1].split(" #")[0]
		elif line[0] == "pop_infile_name":
			pop_infile_name = line[1].split(" #")[0]
		elif line[0] == "ref_iso":
			ref_iso = line[1].split(" #")[0]
		elif line[0] == "ref_contig":
			ref_contig = line[1].split(" #")[0]
		elif line[0] == "len_block_threshold":
			len_block_threshold = int(line[1].split(" #")[0])
		elif line[0] == "gap_prop_thresh":
			gap_prop_thresh = float(line[1].split(" #")[0])

parameter_file.close()

#these directories will generate if they do not already exist
input_dir = project_dir+"align/"
alignment_dir = input_dir+"alignment_blocks/"

#these are output file names
strain_list_filename = "strain_names.txt"
MSA_name = output_prefix+".core.fasta"
phy_prefix = output_prefix
block_loc_filename = output_prefix+".block_location.txt"



#################   MAIN   #################
## Unwrap contig fasta files, also remove all '-' and '.' characters from filenames and sequence headers. MUGSY output will not be readable by these scripts otherwise.
input_file_list = [f for f in os.listdir(input_contig_dir) if f.endswith(contig_extension)]
strain_list_file = open(project_dir+strain_list_filename,"w")
head = ''
for strain_line in input_file_list:
	strain = strain_line.strip().split(contig_extension)[0]
	infile = open(input_contig_dir+strain_line,"r")
	strain = strain.replace(".","").replace("-","")
	strain_list_file.write(strain+"\n")
	outfile = open(contig_dir+strain+".fa","w")
	first_line = True
	contig_no = 0 #counter for contig number (sequential numbering relative to headers in original fasta file)
	for line in infile:
		line = line.strip()
		if line[0] == ">":
			contig_no += 1
			line = line.split(" ")[0]
			head = ">"+strain+"_"+str(contig_no)
			if first_line == True: #just the first line
				outfile.write(head +"\n")
				first_line = False
			else: #every other line
				outfile.write("\n"+ head +"\n")
		else:
			outfile.write(line)
	outfile.write("\n")
infile.close()
outfile.close()
strain_list_file.close()

## Generate MUGSY multiple genome alignment
if os.path.isdir(input_dir) == False:
	os.makedirs(input_dir)

infile = open(project_dir+strain_list_filename,"r")
lis = []
for line in infile:
	line = line.strip()
	line = line + ".fa"
	#line = line.replace("-","")
	lis.append(line)
infile.close()

print(output_prefix)

mugsyline = "mugsy "
for strain in lis:
	mugsyline = mugsyline + contig_dir+strain +" "
mugsyline = mugsyline + "--directory "+input_dir+" --prefix "+output_prefix
print(mugsyline)
#os.system(MUGSY_source+"\n"+mugsyline)
os.system(mugsyline)
