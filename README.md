# bnBreastCancer
This is a program designed to generate projections as to the existence and recurrence of breast cancer

¿En qué consiste?
==============
Es un programa desarrollado en Matlab que utiliza la librería Bayes Net Toolbox para realizar inferencias sobre la red diseñada para el diagnóstico y detección de la recurrencia de cáncer de mama. Como componente extra se diseña una interfaz gráfica con el objetivo de facilitar su uso al usuario.

Características
==============
-	Fácil uso

Requerimientos
==============
El programa está desarrollado en Matlab; requiere de las librerías xml2struct para extraer, del XML generado por Genie, la información necesaria para ser utilizada por la librería Bayes Net Toolbox para realizar las inferencias con la función junction tree como Inferencia exacta para BNs. 

1.	Xml2struct

XML2STRUCT es un programa de MATLAB que lee un archivo que contiene datos XML y produce una estructura MATLAB correspondiente. Desarrollado por Wouter Falkena, ASTI, TUDelft.

2.	Bayes Net Toolbox

BNT admite muchos tipos de distribuciones de probabilidad condicional (nodos)
- Tabular (multinomial)
- Gaussian
- Softmax (logistic/ sigmoid)
- Multi-layer perceptron (neural network)
- Noisy-or
- Deterministic

BNT soporta muchos algoritmos de inferencia diferentes
- Inferencia exacta para BNs estáticos:
-	Inferencia aproximada para BNs estáticos:
-	Inferencia exacta para DBNs:
o	junction tree
-	Inferencia aproximada para DBNs:

BNT soporta varios métodos para el aprendizaje de parámetros, BNT soporta varios métodos para el aprendizaje de estructuras
Para más información sobre el manejo de la instalación, creación de la red bayes y generación de inferencias con la librería bnt, así como comprender la funcionalidad de la librería xml2struct, se dejan los siguientes enlaces.

Enlace [xml2strct](https://github.com/azag0/xml2struct)

Enlace [bnt](http://bayesnet.github.io/bnt/docs/usage.html#basics)

Instalación de componentes y ejecución del programa
==============
1.	Descargar la librería Bayes Net Toolbox desde el repositorio github: https://github.com/bayesnet/bnt y la librería xml2struct desde el repositorio github: https://github.com/azag0/xml2struct.
2.	Ubicar el SetPath de Matlab para incorporar las librerías y estas puedan ser utilizada en cualquier lugar en donde exista un archivo de Matlab en el sistema.
3.	Ubicar en el sistema donde se encuentra el programa desarrollado “bntProyectoV2.m” y “bntProyectoV2.fig”, junto a la red exportada de Genie con el nombre “TRAINED.xdsl”.
4.	Ejecutar en la ventana de comando de window de Matlab bntProyectoV2.

¿Cómo funciona?
==============
Una vez se está situado en la pantalla principal del lado de las Ocurrencias, se generan las probabilidades de la forma P(X) y al adicionar evidencias se generan las probabilidades de la forma P(X|Y).

