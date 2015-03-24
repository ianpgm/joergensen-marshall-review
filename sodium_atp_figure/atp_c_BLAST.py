import os
from Bio import SeqIO

#BLAST Peru Margin ATP synthase c subunit sequences against NR database

faa_query_sequences = open("real_PM_c_subunits.faa")

blast_result_dict = dict()

output_file = open("real_PM_c_subunits_vs_NR.tdt","w")

for record in SeqIO.parse(faa_query_sequences, "fasta"):
	single_sequence_query = open("temp.faa", "w")
	SeqIO.write(record, single_sequence_query, "fasta")
	os.system("./blastp -remote -db nr -max_target_seqs 1 -outfmt 7 -query temp.faa -out temp_results.tdt")
	for line in open("temp_results.tdt"):
		if line.startswith("#"):
			continue
		else:
			output_file.write(line)



