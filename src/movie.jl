# movie.jl

using JSON3

contenido_fichero = read("src/posiciones.json", String)
datos_recuperados = JSON3.read(contenido_fichero)
x1 = datos_recuperados["x1"]
x2 = datos_recuperados["x2"]

