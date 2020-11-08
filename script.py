from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

# Definición de la función summarize_contents
def summarize_contents(filename):
	listaRuta = os.path.split(filename)
	extension = os.path.splitext(filename)
	# Tipo de archivo genbank o fasta 
	if(extension[1] == ".gbk"):
		tfile = "genbank"
	else: 
		tfile = "fasta"
	records = list(SeqIO.parse(filename, tfile))
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
	for seq_record in SeqIO.parse(filename, tfile):
		diccionario['names:'].append(seq_record.name)
		diccionario['IDs:'].append(seq_record.id)
		diccionario['descriptions:'].append(seq_record.description)
	return diccionario

# Definición de la función concatenate_and_get_reverse_of_complement()
def concatenate_and_get_reverse_of_complement(sec1, sec2):
	# Concatena las 2 secuencias y las convierte en objeto Seq
	concatenado = Seq(sec1 + sec2)
	# Devuelve el reverso complementario en mayúsculas de la concatenación
	return concatenado.reverse_complement().upper()

# Llamada a la función summarize_contents
if __name__ == "__main__":
	filename =  os.path.abspath("data/ls_orchid.gbk")
	secuencia_1 = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
	secuencia_2 = "AGTCAG"
	resultado = summarize_contents(filename)
	print("Resultado summarize_contents:\n", resultado)
	resultado = concatenate_and_get_reverse_of_complement(secuencia_1, secuencia_2)
	print("\nResultado concatenate_and_get_reverse_of_complement:\n", resultado)