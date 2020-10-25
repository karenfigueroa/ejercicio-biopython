from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

# Archivo .gbk de mi escritorio, cambiar a la dirección del nuevo archivo .gbk a leer
filename =  "/mnt/c/Users/karen/Desktop/BIOINFORMATICA/biopython-notebook/notebooks/data/ls_orchid.gbk"

# Definición de la función summarize_contents
def summarize_contents(filename):
	listaRuta = []
	listaRuta = os.path.split(filename)
	records = list(SeqIO.parse(filename, "genbank"))
	# Crear diccionario 
	diccionario = {}
	# Nombre del archivo, ruta y número de registros
	diccionario['file:'] = listaRuta[1]
	diccionario['path:'] = listaRuta[0]
	diccionario['num_records:'] = len(records)
	# Crear listas en diccionario
	diccionario['names:'] = []
	diccionario['IDs:'] = []
	diccionario['descriptions:'] = []
	# Registros
	for seq_record in SeqIO.parse(filename, "genbank"):
		diccionario['names:'].append(seq_record.name)
		diccionario['IDs:'].append(seq_record.id)
		diccionario['descriptions:'].append(seq_record.description)
	return diccionario

# Llamada a la función summarize_contents
if __name__ == "__main__":
	resultado = summarize_contents(filename)
	print(resultado)