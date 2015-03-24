import os
import urllib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord

#Function to remove gaps from an record contained an aligned sequence (for the output of "real" ATP synthase c subunits, i.e. those that don't have a gap at the Ilyobacter tartaricus position 67)
def remove_gaps(record):
	new_seq = str()
	for char in record.seq:
		if char == "-":
			continue
		else:
			new_seq = new_seq + char
	new_record = record
	new_record = SeqRecord(Seq(new_seq, generic_protein), id=record.id, description = record.description)
	return(new_record)

#Read in the resulting alignment, determine the status of each sequence based on I. tartaricus pos. 67
aligned_c_subunits = open("PM_c_subunits_IT_aligned.stockholm")

sodium_residues = ['T','S','R','Q']
proton_residues = ['L','I','V','A','F']
ion_dict = dict()

real_c_subunits = list()
misassigned_subunits = list()

for record in SeqIO.parse(aligned_c_subunits, "stockholm"):
	IT_pos67 = record.seq[365]
	print(record.id)
	print(IT_pos67)
	if IT_pos67 in sodium_residues:
		ion_dict[record.id] = "Na"
		real_c_subunits.append(remove_gaps(record))
		continue
	if IT_pos67 in proton_residues:
		ion_dict[record.id] = "H"
		real_c_subunits.append(remove_gaps(record))
		continue
	if IT_pos67 == "-":
		ion_dict[record.id] = "non-c-subunit"
		misassigned_subunits.append(remove_gaps(record))
		continue
	else:
		ion_dict[record.id] = "other"
		real_c_subunits.append(remove_gaps(record))
		print(IT_pos67)

#Output of 'real' ATP synthase c subunits (i.e. those I don't think have been misassigned by HMMER)
output_fasta = open("real_PM_c_subunits.faa", "w")
SeqIO.write(real_c_subunits, output_fasta, "fasta")

#Output of misassigned ATP synthase c subunits (i.e. those I think have been misassigned by HMMER)
mis_output_fasta = open("misassigned_PM_c_subunits.faa", "w")
SeqIO.write(misassigned_subunits, mis_output_fasta, "fasta")


#Calculate the abundance of each type in each sample based on coverage in each sequence ID
result_dict = dict()

for sequence in ion_dict:
	if sequence == "AF522463_3":
		continue
	sample = sequence.split("_")[0]
	coverage = int(sequence.split("[")[1].split("]")[0].split("=")[1])
	atp_synt_type = ion_dict[sequence]
	if sample not in result_dict.keys():
		result_dict[sample] = dict()
	if atp_synt_type not in result_dict[sample].keys():
		result_dict[sample][atp_synt_type] = int()
	result_dict[sample][atp_synt_type] = result_dict[sample][atp_synt_type] + coverage

print(result_dict)

print('sample\tH+\tNa+\tother\tnon-c-subunit')

headers = ['Na','H','other','non-c-subunit']

for sample in result_dict:
	for header in headers:
		if header not in result_dict[sample].keys():
			result_dict[sample][header] = 0	
	print(sample + "\t" + str(result_dict[sample]['H']) + "\t" + str(result_dict[sample]['Na']) + "\t" +str(result_dict[sample]['other']) + "\t" +str(result_dict[sample]['non-c-subunit']))
