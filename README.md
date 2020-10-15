ejercicio-biopython

Descripción de script:
Contiene la importación de las bibliotecas de Biopython y una función llamada summarize_contents que imprime un resumen del contenido de un archivo genbank (.gbk) con la siguiente estructura:
file: [nombre de archivo]
path: [ruta al archivo]
num_records: [numero de registros]
records:
		- id: [id del primer registro]
		  name: [nombre]
		  description: [descripción]
		- id: [id del segundo registro]
		  name: [nombre]
		  description: [descripción]
		- ...

Descripción función summarize_contents:
La función recibe el nombre de un archivo que se asume que está en formato genbank y con os.path.split, se utiliza una lista (listaFile) que guarda la ruta y el nombre del archivo por separado para su posterior impresión. Después se imprime el número de records, y cada uno de ellos con un ciclo for que itera hasta el último registro. Imprimiendo su id, nombre y descripción. 
