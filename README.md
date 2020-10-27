ejercicio-biopython

Descripción de script:
Contiene la importación de las bibliotecas de Biopython, una función llamadada summarize_contents que devuelve como resultado un resumen del contenido de un archivo genbank (.gbk) o fasta (.fasta), el cual es almacenado en una estructura de  directorio de python. Por último, se imprime el directorio resultante. 

Descripción función summarize_contents:
La función recibe el nombre de un archivo y con os.path.split, se utiliza una lista (listaFile) que guarda la ruta y el nombre del archivo por separado. Con os.path.splitext, obtenemos el tipo de archivo que se esta leyendo, genbank o fasta, y se guarda en extension. De acuerdo al valor de esta variable, se asigna a la variable tfile, el tipo de archivo que se esta leyendo. Posteriormente, se crea un diccionario, donde se guarda cada característica del archivo leído. La función summarize_contents devuelve un directorio conteniendo el file, path, num_records y cada records con su respectivo id, nombre y descripción.

Descripción función test_script:
Contiene la importación de las bibliotecas utilizadas para llevar a cabo la prueba de unidad de la función summarize_contents con la función test_summarize_contents, contenida dentro de la clase MiPrueba

Descripción función summarize_contents:
La función prueba que dos diccionarios son iguales respectivamente, primero guarda el output de cada archivo genbank y fasta contenido en la carpeta data, después guarda en otra variable el resultado de la llamada a la función summarize_contents con os.path.abspath("data/ls_orchid.gbk") como parámetro de la función. Por último, asegura que sean iguales con assertDictEqual. 
