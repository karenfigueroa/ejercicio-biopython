ejercicio-biopython

Descripción de script:
Contiene la importación de las bibliotecas de Biopython y una función llamada summarize_contents que imprime un resumen del contenido de un archivo genbank (.gbk)

Descripción función summarize_contents:
La función recibe el nombre de un archivo que se asume que está en formato genbank, y utilizando librerias de biopython, se imprime un resumen del contenido del archivo con la siguiente estructura:

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
