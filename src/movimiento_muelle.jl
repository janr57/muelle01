# movimiento_muelle.jl

using JSON3
using GLMakie

# 1) Recuperación de los arrays y convertirlos en tuplas
contenido_fichero = read("src/posiciones.json", String)
datos_recuperados = JSON3.read(contenido_fichero)
x1_array = datos_recuperados["x1"]
x2_array = datos_recuperados["x2"]
x1_tuple = Tuple(x1_array)
x2_tuple = Tuple(x2_array)

# Valores máximos y mínimos
x1min = minimum(x1)
x1max = maximum(x1)
x2min = minimum(x2)
x2max = maximum(x2)

xmin = floor(min(x1min, x2min))
xmax = 2 * ceil(max(x1max, x2max))
ymin = -0.02
y = 0.0
ymax = 0.101

# 2) Creación de la figura y los ejes
f = Figure()
ax = Axis(f[1,1], title = "Muelle en movimiento")

# Límites
limits!(ax, xmin - 0.015 * (xmax - xmin), xmax + 0.015 * (xmax - xmin), ymin, ymax)

# 2)  Definición de las posiciones x1 y x2 como observables
x1 = x1_tuple[1]
x2 = x2_tuple[1]
pos1 = Observable(Point2f(x1, y))
pos2 = Observable(Point2f(x2, y))

# 3) Representación de las masas
lines!(ax, [x1, x2], [y, y], color = :black)
lines!(ax, pos2, color = :black)
scatter!(ax, pos1, markersize = 20, color = :red)
scatter!(ax, pos2, markersize = 20, color = :red)
f




