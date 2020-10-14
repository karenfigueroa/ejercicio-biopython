from Bio.Seq import Seq 
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import os

def summarize_contents(filename):
	record = SeqIO.read(filename, "genbank")
	print("Name: ", record.name)
	print("Path: ", os.path.dirname(filename))
	records = list(SeqIO.parse (filename, "genbank"))
	print("num_records = %i records" % len(records))
	
	for seq_record in SeqIO.parse(filename, "genbank"):
		print("ID: ", record.id)
		#Pendiente location
	
summarize_contents(filename)
