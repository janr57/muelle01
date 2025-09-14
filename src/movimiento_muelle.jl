# movimiento_muelle.jl

using JSON3
using GLMakie

# 1) Recuperación de los arrays y convertirlos en tuplas
contenido_fichero = read("src/posiciones.json", String)
datos_recuperados = JSON3.read(contenido_fichero)
# Estos valores recuperados están normalizados de 0 a 100
x1_array = datos_recuperados["x1"]
x2_array = datos_recuperados["x2"]
x1_tuple = Tuple(x1_array)
x2_tuple = Tuple(x2_array)

len = length(x1_tuple)

# 2) Creación de la figura y los ejes
# Límites de los ejes
incx = 1.6
xmin = 0.0 - incx
xmax = 100.0 + incx
ymin = -0.02
ymax = 0.101

# Coordenada y de las masas en la figura
y = 0.0

f = Figure()
ax = Axis(f[1,1], title = "MUELLE EN MOVIMIENTO")

# Límites
limits!(ax, xmin, xmax, ymin, ymax)

# 2)  Definición de las posiciones x1 y x2 como observables
x1_ini = x1_tuple[1]
x2_ini = x2_tuple[1]
pos1 = Observable(Point2f(x1_ini, y))
pos2 = Observable(Point2f(x2_ini, y))
muelle = Observable((pos1, pos2))

# 3) Representación de las masas
#lines!(ax, [x1_ini, x2_ini], [y, y], color = :red)
lines!(ax, [muelle[][1][][1], muelle[][2][][1]], [0.0, 0.0], color = :black)
scatter!(ax, pos1, markersize = 20, color = :black)
scatter!(ax, pos2, markersize = 20, color = :black)

display(f)

# 4) Animación del sistema
for i in 2:len
    # Nuevas coordenadas
    x1_new = x1_tuple[i]
    x2_new = x2_tuple[i]
    
    # Nuevos valores de los observables
    pos1[] = Point2f(x1_new, y)
    pos2[] = Point2f(x2_new, y)
    muelle[] = (pos1, pos2)
    #notify!(muelle)

    #lines!(ax, [x1_new, x2_new], [y, y], linewidth = 1, color = :blue)

    sleep(0.005)
end

println("FIN")



