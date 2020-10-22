from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# Archivo .gbk de mi escritorio, cambiar a la dirección del nuevo archivo .gbk a leer
filename =  "/mnt/c/Users/karen/Desktop/BIOINFORMATICA/biopython-notebook/notebooks/data/ls_orchid.gbk"

# Definición de la función summarize_contents
def summarize_contents(filename):
	listaRuta = []
	listaRuta = os.path.split(filename)
	cadena = " "
	# File y ruta
	cadena = ("file: "+ listaRuta[1] + "\npath: " + listaRuta[0])
	# Número de registros
	all_records=[]
	records = list(SeqIO.parse(filename, "genbank"))
	cadena += ("\nnum_records: " + str(len(records)))
	cadena += ("\nrecords:")
	# Registros
	for seq_record in SeqIO.parse(filename, "genbank"):
		all_records.append(seq_record.name)
		cadena += ("\n- id:" + str(seq_record.id))
		cadena += ("\nname: " + seq_record.name)
		cadena += ("\ndescription: " + str(seq_record.description))
		cadena += ("\n")
	return cadena
# Llamada a la función summarize_contents
resultado = summarize_contents(filename)
print(resultado)
