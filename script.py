from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# Archivo .gbk de mi escritorio, cambiar a la dirección del nuevo archivo .gbk a leer
filename =  "/mnt/c/Users/karen/Desktop/BIOINFORMATICA/biopython-notebook/notebooks/data/ls_orchid.gbk"

# Definición de la función summarize_contents
def summarize_contents(filename):
        listaRuta = []
        listaRuta = os.path.split(filename)
        # File y ruta
        print("file:", listaRuta[1], "\npath:", listaRuta[0])
        #print ("path: ", os.path.dirname(filename))
        all_records=[]
        records = list(SeqIO.parse(filename, "genbank"))
        # Número de registros
        print("num_records = %i records" % len(records))
        print("records:")
        # Registros
        for seq_record in SeqIO.parse(filename, "genbank"):
                all_records.append(seq_record.name)
                print("- id:",seq_record.id)
                print("name: ", seq_record.name)
                print("description: ", seq_record.description)
                print("\n")
        
# Llamada a la función summarize_contents
summarize_contents(filename)
