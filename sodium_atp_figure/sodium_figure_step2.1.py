import os
import urllib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord

sequence_names = set()

for line in open("hmmscan_output_table"):
	if line.startswith("#"):
		continue
	else:
		sequence_names.add(line[32:].split(" ")[0])

#Retrieve sequences of identified c subunit genes
PM_c_subunits = open("PM_c_subunits.faa", "w")

for record in SeqIO.parse(open("Peru_Margin_proteins.faa"), "fasta"):
	if record.id in sequence_names:
		PM_c_subunits.write(">" + record.description + "\n" + str(record.seq) + "\n")
		

#Add the Ilyobacter tartaricus c subunit gene to the file for alignment
os.system("cat PM_c_subunits.faa I_tartaricus_c_subunit.fasta > PM_c_subunits_IT.faa;hmmalign -o PM_c_subunits_IT_aligned.stockholm ATP-synt_C.hmm PM_c_subunits_IT.faa")
