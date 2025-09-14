# movimiento_muelle.jl

using JSON3
using GLMakie

# 1) Recuperación de los arrays y convertirlos en tuplas
contenido_fichero = read("src/posiciones.json", String)
datos_recuperados = JSON3.read(contenido_fichero)
# Estos valores recuperados están normalizados de 0 a 100
xp1_array = datos_recuperados["xp1"]
xp2_array = datos_recuperados["xp2"]
xp1_tuple = Tuple(xp1_array)
xp2_tuple = Tuple(xp2_array)

len = length(xp1_tuple)

# 2) Creación de la figura y los ejes
# Límites de los ejes para poder representar las masas de forma completa
# Margen 
incxp = 1.1
xpmin = 0.0 - incxp
xpmax = 100.0 + incxp
ypmin = 0.0 - 0.19
ypmax = 0.0 + 1.6

# Coordenada y de las masas en la figura
yp = 0.0

f = Figure(size = (800, 150))
ax = Axis(f[1,1], title = "MUELLE EN MOVIMIENTO")

# Límites
limits!(ax, xpmin, xpmax, ypmin, ypmax)

index = Observable(1)
pos1 = @lift(Point2f(xp1_tuple[$index], yp))
pos2 = @lift(Point2f(xp2_tuple[$index], yp))
#muelle = @lift((Point2f(xp1_tuple[$index], yp), Point2f(xp2_tuple[$index], yp)))
muelle = @lift([xp1_tuple[$index], xp2_tuple[$index]])

# Representa el muelle mediante una línea negra
#lines!(ax, [pos1[][1], pos2[][1]], [yp, yp], linewidth = 1, color = :black)
# Representa las masas 1 y 2
lines!(ax, muelle, [yp, yp], linewidth = 1, color = :red)
scatter!(ax, pos1, markersize = 20, color = :black)
scatter!(ax, pos2, markersize = 20, color = :black)

# Anima la masa 1
# Muestra la gráfica inicial
display(f)

for i in 2:len
    # Calculate new coordinates
    # Actualiza la posición de la masa 1 con las nuevas coordenadas
    pos1[] = Point2f(xp1_tuple[i], yp)
    pos2[] = Point2f(xp2_tuple[i], yp)
    #muelle[] = (Point2f(xp1_tuple[i], yp), Point2f(xp2_tuple[i], yp))
    muelle[] = [xp1_tuple[i], xp2_tuple[i]]

    sleep(0.001)
end

println("FIN")



