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
	diccionario['mRNA'] = secuencia.transcribe().upper()
	diccionario['proteins'] = []
	diccionario['stop_codons'] = []
	
	# Proceso de búsqueda de proteínas y stop codons
	aminoacidos = secuencia.translate(table = num_table)	
	# Exportando codones de inicio y de parada de la tabla estandar
	codons_table = CodonTable.unambiguous_dna_by_id[num_table]
	start_codons_list = codons_table.start_codons
	stop_codons_list = codons_table.stop_codons
	start_codon, stop_codon = False, False
	i = 0
	while i < len(aminoacidos):
		if secuencia[i*3:i*3+3].upper() in start_codons_list:
			start_codon = True
			if i+1 == len(aminoacidos):
				diccionario['proteins'].append(aminoacidos[i:])
				break
			j = i+1
			while j < len(aminoacidos):
				if secuencia[j*3:j*3+3].upper() in stop_codons_list:
					stop_codon = True
					diccionario['proteins'].append(aminoacidos[i:j])
					diccionario['stop_codons'].append(secuencia[j*3:j*3+3].upper())
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

# Definición de la función print_proteins_and_codons_using_mitocondrial_yeast_table()
def print_proteins_and_codons_using_mitocondrial_yeast_table(sec):
	# Convertir secuencia en objeto Seq
	secuencia = Seq(sec)	
	# Crea diccionario
	diccionario = {}
	# Guarda el diccionario resuelto por la función solve_dictionary
	diccionario = solve_dictionary(secuencia, 3)
	return diccionario

#Definición de la función extract_sequences():
def extract_sequences(file, format):
	# Parte 2 ejercicio 4
	#direccion = os.path.abspath(file)
	#records = list(SeqIO.parse(direccion, "fasta"))
	#for i in range(len(records)):
	#	filename = open(f"sequence{i+1}.fasta", "w")
	#	filename.write('>' + records[i].id + os.linesep)
	#	filename.write(str(records[i].seq))
	#	filename.close()

	# Parte 3 ejercicio 4
	root_ext = os.path.splitext(file)
	if root_ext[1] == ".fasta" and format.lower() == "genbank":  
		# Convierte archivo fasta a genbank
		SeqIO.convert(file, "fasta", "my_example.gbk", "genbank", molecule_type = "DNA")
		direccion = os.path.abspath("my_example.gbk")
		records = list(SeqIO.parse(direccion, "genbank"))
		for i in range(len(records)):
			filename = open(f"sequence{i+1}.gbk", "w")
			filename.write(str(records[i].format("genbank")))
			filename.close()
		os.remove("my_example.gbk")
	else: 
		# Lanza error ya que sólo se reciben archivos .fasta y se generan en formato .gbk
		return "Error: this program can only convert .fasta files to N files in .genbank format"

# Definición de la función extract_sequences_revcomp():
def extract_sequences_revcomp(file):
	root_ext = os.path.splitext(file)
	if root_ext[1] == ".fasta":
		direccion = os.path.abspath(file)
		records = list(SeqIO.parse(direccion, "fasta"))
		filename = open("reverse_complement.fasta", "w")
		for i in range(len(records)):
			filename.write(str(records[i].reverse_complement(id = True).format("fasta")))
		filename.close()
	else:
		return "Error: the only file format allowed is .fasta"

# Llamada a las funciones
if __name__ == "__main__":
	filename =  os.path.abspath("data/ls_orchid.gbk")
	secuencia_1 = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
	secuencia_2 = "ATAATAATGGTGGGGTAATAGTAAGTGAAAAAATAGCTTTTT"
	resultado = summarize_contents(filename)
	print("Resultado summarize_contents:\n", resultado)

	resultado = concatenate_and_get_reverse_of_complement(secuencia_1, secuencia_2)
	print("\nResultado concatenate_and_get_reverse_of_complement:\n", resultado)

	resultado = print_protein_and_codons_using_standard_table(secuencia_2)
	print("\nResultado print_protein_and_codons_using_standard_table:\n", resultado)
	
	resultado = print_proteins_and_codons_using_mitocondrial_yeast_table(secuencia_2)
	print("\nResultado print_proteins_and_codons_using_mitocondrial_yeast_table:\n", resultado)

	extract_sequences("data/sequences.fasta", "genbank")

	extract_sequences_revcomp("data/sequences.fasta")