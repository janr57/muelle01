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
incxp = 1.6
xpmin = 0.0 - incx
xpmax = 100.0 + incx
ypmin = -0.02
ypmax = 0.101

# Coordenada y de las masas en la figura
yp = 0.0

f = Figure()
ax = Axis(f[1,1], title = "MUELLE EN MOVIMIENTO")

# Límites
limits!(ax, xpmin, xpmax, ypmin, ypmax)

index = Observable(1)
pos1 = @lift(Point2f(xp1_tuple[$index], yp))
pos2 = @lift(Point2f(xp2_tuple[$index], yp))

# Representa la masa 1
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

    sleep(0.001)
end








#pos2 = @lift(Point2f(x2_tuple[$index], y))

# Step 2: Use @lift to define the point's position
#point_pos = @lift(Point2f(cos(deg2rad($angle)), sin(deg2rad($angle))))

# Step 3: Plot the point
#scatter!(ax, point_pos, markersize = 20, color = :red)






# 2)  Definición de las posiciones x1 y x2 como observables
#x1_ini = x1_tuple[1]
#x2_ini = x2_tuple[1]



#pos1 = Observable(Point2f(x1_ini, y))
#pos2 = Observable(Point2f(x2_ini, y))
#muelle = Observable((pos1, pos2))


## Step 2: Use @lift to define the point's position
#muelle = @lift(Point2f(cos(deg2rad($angle)), sin(deg2rad($angle))))

## 3) Representación de las masas
##lines!(ax, [x1_ini, x2_ini], [y, y], color = :red)
#lines!(ax, [muelle[][1][][1], muelle[][2][][1]], [0.0, 0.0], color = :black)
#scatter!(ax, pos1, markersize = 20, color = :black)
#scatter!(ax, pos2, markersize = 20, color = :black)
#
#display(f)
#
## 4) Animación del sistema
#for i in 2:len
#    # Nuevas coordenadas
#    x1_new = x1_tuple[i]
#    x2_new = x2_tuple[i]
#    
#    # Nuevos valores de los observables
#    pos1[] = Point2f(x1_new, y)
#    pos2[] = Point2f(x2_new, y)
#    muelle[] = (pos1, pos2)
#    #notify!(muelle)
#
#    #lines!(ax, [x1_new, x2_new], [y, y], linewidth = 1, color = :blue)
#
#    sleep(0.005)
#end

println("FIN")



