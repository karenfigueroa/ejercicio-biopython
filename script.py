from Bio.Seq import Seq 
from Bio.SeqFeatures import SeqFeatures, FeaturesLocation
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
def summarize_contents(filename):
	record = SeqIO.read("filename", "genbank")
	print(record)

