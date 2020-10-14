from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import os

# Archivo .gbk de mi escritorio, cambiar a la dirección del nuevo archivo .gbk a leer
filename =  "/mnt/c/Users/karen/Desktop/BIOINFORMATICA/biopython-notebook/notebooks/data/ls_orchid.gbk"

# Definición de la función summarize_contents
def summarize_contents(filename):
        all_records=[]
        records = list(SeqIO.parse(filename, "genbank"))
        print ("Path: ", os.path.dirname(filename))
        print("num_records = %i records" % len(records))
        print("\n\n")
        for seq_record in SeqIO.parse(filename, "genbank"):
                all_records.append(seq_record.name)
                print("Name: ", seq_record.name)
                print("ID:",seq_record.id)
                print("Location:")
                for seq_feature in seq_record.features:
                        print('Start: %d, Stop: %d'%(int(seq_feature.location.start),int(seq_feature.location.end)))
                print("\n")
        
# Llamada a la función summarize_contents
summarize_contents(filename)
