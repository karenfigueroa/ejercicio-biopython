from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
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

#Definición de la función que llena los datos del diccionario con mRNA, proteínas y codones de paro de una secuencia
def solve_dictionary(secuencia, num_table):
	# Crea diccionario
	diccionario = {}
	# mRNA, proteínas y codones de paro
	diccionario['mRNA'] = secuencia.transcribe()
	diccionario['proteins'] = []
	diccionario['stop_codons'] = []
	
	# Proceso de búsqueda de proteínas y stop codons
	aminoacidos = secuencia.translate(table = num_table)	
	# Exportando codones de inicio y de parada de la tabla estandar
	codons_table = CodonTable.unambiguous_dna_by_id[num_table]
	start_codons_list = codons_table.start_codons
	stop_codons_list = codons_table.stop_codons
	start_codon, stop_codon = False, False
	pos_inicio, pos_final = None, None
	i = 0
	while i < len(aminoacidos):
		if secuencia[i*3:i*3+3].upper() in start_codons_list:
			start_codon = True
			pos_inicio = i
			if i+1 == len(aminoacidos):
				diccionario['proteins'].append(aminoacidos[i:])
				break
			j = i+1
			while j < len(aminoacidos):
				if secuencia[j*3:j*3+3].upper() in stop_codons_list:
					stop_codon = True
					pos_final = j
					diccionario['proteins'].append(aminoacidos[i:j])
					diccionario['stop_codons'].append(secuencia[j*3:j*3+3])
					start_codon, stop_codon = False, False
					i = j
					break
				j += 1
		if(start_codon == True and stop_codon == False):
			diccionario['proteins'].append(aminoacidos[i:])
			break
		i += 1

	# Verifica si no se encontraron proteínas en la secuencia dada
	if (diccionario['proteins'] == []):
		diccionario['proteins'] = "Not found proteins"
	if diccionario['stop_codons'] == []:
		diccionario['stop_codons'] = "Not found stop codons"
	
	return diccionario

# Definción de la función print_protein_and_codons_using_standard_table()
def print_protein_and_codons_using_standard_table(sec):
	# Convertir secuencia en objeto Seq
	secuencia = Seq(sec)	
	# Crea diccionario
	diccionario = {}
	# Guarda el diccionario resuelto por la función solve_dictionary
	diccionario = solve_dictionary(secuencia, 1)
	return diccionario

# Llamada a las funciones
if __name__ == "__main__":
	filename =  os.path.abspath("data/ls_orchid.gbk")
	secuencia_1 = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
	secuencia_2 = "ATG"
	resultado = summarize_contents(filename)
	print("Resultado summarize_contents:\n", resultado)

	resultado = concatenate_and_get_reverse_of_complement(secuencia_1, secuencia_2)
	print("\nResultado concatenate_and_get_reverse_of_complement:\n", resultado)

	resultado = print_protein_and_codons_using_standard_table(secuencia_1)
	print("\nResultado print_protein_and_codons_using_standard_table:\n", resultado)
