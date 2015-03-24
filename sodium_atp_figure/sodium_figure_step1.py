import os
import urllib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.SeqRecord import SeqRecord

#Retrieve Orsi et al. 2013 Peru Margin metatranscriptomic data
base_url = "http://metagenomics.anl.gov/metagenomics.cgi?page=DownloadMetagenome&stage=350&action=download&file=350.genecalling.coding.faa.gz&metagenome=" 
MG_RAST_accessions = ['4515478.3', '4515477.3', '4515476.3', '4510337.3', '4510336.3', '4510335.3']

for accession in MG_RAST_accessions:
	urllib.urlretrieve(base_url + accession, accession + "_PM.faa.gz")
	os.system("gunzip " + accession + "_PM.faa.gz")
	
os.system("cat *_PM.faa > Peru_Margin_proteins.faa")
	
#Retrieve PFAM hidden markov model (HMM) descrbing ATP synthase c subunit
urllib.urlretrieve("http://pfam.sanger.ac.uk/family/PF00137/hmm","ATP-synt_C.hmm")

#Make binary hidden markov model for ATP synthase c subunits
os.system("hmmpress ATP-synt_C.hmm")

#Scan for genes encoding ATP synthase c subunit
os.system("hmmscan --tblout hmmscan_output_table ATP-synt_C.hmm Peru_Margin_proteins.faa")