# muelle01

Objetivo:
--------
Analizar el movimiento de dos masas *m1* y *m2*, unidas a un muelle de 
constante *k*, cuando se encuentra comprimido con una cierta elongación *x0* y la
masa *m1* se encuentra contra una pared vertical, desde el momento en que se 
suelta. Se desprecia el rozamiento con el suelo.

Procedimiento:
-------------
Se obtienen las ecuaciones diferenciales de las dos fases del movimiento del sistema:
1.	Mientras la masa *m1* se encuentra en reposo contra la pared y solo se
	desplaza *m2*.

2. Cuando *m1* se libera de la pared y comienza a oscilar junto con *m2*.

Se resuelven analítica y numéricamente las ecuaciones y se representan
gráficamente, en función del tiempo, para cada una de las masas y del centro 
de masas la:
- Posición.
- Velocidad.
- Aceleración.
- Energía cinética, potencial y mecánica.

Además se representa la elongación del muelle, también en función del tiempo.

Lenguaje y módulos:
-------------------
Se utiliza el lenguaje Julia, junto con los paquetes Pluto.jl, PlutoUI.jl,
DifferentialEquations.jl y CairoMakie.jl

Puesta en marcha:
------------------

- SI ES LA PRIMERA VEZ QUE SE EJECUTA EL PROGRAMA EN ESE ORDENADOR:
	1. Comprobar que tenemos la última versión de *Julia*:
		$ juliaup update
	2. Entrar en Julia REPL
		$ julia
	3. Instalar Pluto
		julia> import Pkg
		julia> Pkg.add("Pluto")
	4. Salir de Julia REPL pulsando CONTROL-d
	5. Entrar en el directorio y activar el proyecto:
		$ cd ~/.../muelle01
		$ julia --project
	6. Entrar en la gestión de paquetes e instanciar el proyecto y salir de la
	   gestión de paquetes (puede tardar un poco):
		julia> ]
		(muelle01) pkg> instantiate
		(muelle01) pkg> (pulsa la tecla DEL, encima de INTRO)
	7. Arrancar Pluto.jl
		julia> import Pluto; Pluto.run()
	8. Abrir el cuaderno de Pluto *src/muelle01.jl* en la sección "Open a notebook" 
	del cuaderno y hacer click en "Open".
	9. Pulsar el botón "Run notebook code" de la parte superior del cuaderno.
	En el cuadro "Status" de la parte inferior del cuaderno, se muestra el avance
	de la primera precompilación (que puede tardar un poco la primera vez).
		

- EN CASO CONTRARIO:
	1. Entrar en el directorio y activar el proyecto:
		$ cd ~/.../muelle01
		$ julia --project=.
	2. Arrancar Pluto:
		julia> import Pluto; Pluto.run()
	3. Abrir el cuaderno de Pluto *src/muelle01.jl* en la sección "My work" 
	del cuaderno y hacer click en "muelle01.jl", o bien, abrir el cuaderno de Pluto
	*src/muelle01.jl* en la sección "Open a notebook" del cuaderno y hacer click en
	"Open".
	4. Pulsar el botón "Run notebook code" de la parte superior del cuaderno.
	En el cuadro "Status" de la parte inferior del cuaderno, se muestra el avance
	de la primera precompilación (que puede tardar un poco la primera vez).




