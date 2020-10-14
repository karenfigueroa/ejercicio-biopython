from Bio.Seq import Seq 
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
def summarize_contents(filename):
	record = SeqIO.read(filename, "genbank")
	print(record)
summarize_contents(filename)

