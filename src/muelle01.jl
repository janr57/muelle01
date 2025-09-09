### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ e89a98ac-917b-4094-871d-6f75a344567b
using DifferentialEquations, Plots,PlutoUI

# ╔═╡ e5209b22-863c-11f0-2e47-e157a5f44a52
md"# MUELLE CON DOS MASAS, COMPRIMIDO CONTRA UNA PARED"

# ╔═╡ 29d4fb45-0f1b-424f-8de6-ed277df6931f
md"## Objetivos planteados con este cuaderno"

# ╔═╡ ad8ad00a-0f02-40bd-9329-e0f76c92bb92
md"""
1. Aprender a resolver ecuaciones diferenciales acopladas sencillas, utilizando el lenguaje `Julia` y los paquetes `Pluto.jl`. `PlutoUI.jl`.
2. Representar gráficamente mediante el paquete ```Plots.jl```.
2. Comparar las soluciones anteriores con un posible análisis analítico.
"""

# ╔═╡ 88a02306-6b34-48df-ad4f-0f232b5ec826
md"## Movimiento del muelle"

# ╔═╡ 29fa4220-9626-4a87-bd99-49ac4ed74df4
md"""
- Tiene lugar en una dimensión (eje ``x``), por tanto, podremos representar los vectores correspondientes mediante escalares, prescindiendo del vector unitario correspondiente, sin que se produzca ambigüedad alguna.
- El origen de coordenadas, ``O``, se elige en la pared (a la izquierda).
"""

# ╔═╡ 8df61bb0-70ee-4744-b0b9-8a71ac519d10
md"### Fases del movimiento"

# ╔═╡ bacf9614-5d5c-4d6b-a2b3-3bb4c7dd8e9f
md"""
1. El muelle empuja a la masa ``m₁`` contra la pared, y esta no se mueve. Esta fase termina cuando el muelle adquiere su longitud natural ``L₀`` y la masa ya no presiona la pared. Como consecuencia, la masa ``m₁`` se libera de ella.
2. En esta fase, el muelle oscila y las dos masas, ``m₁`` y ``m₂`` también lo hacen. Además, el centro de masas del sistema se desplaza con velocidad constante.
"""


# ╔═╡ cc7fe116-32e5-45da-b0c7-0fff78039593
md"### Primera fase"

# ╔═╡ 96c1c896-4636-4a37-a6b6-40c4c1c6e729
md"#### Parámetros"

# ╔═╡ b8943312-9df6-4fdd-a94e-31ccff232b63
md"""
- ``m₁ \equiv \text{ Masa situada contra la pared (a la izquierda del sistema) } (kg)``. En esta primera fase no se separa de la pared.
- ``m₂ \equiv \text{ Masa situada a cierta distancia de la pared (a la derecha del sistema) } (kg)``. Esta masa se mueve sujeta a la fuerza que le ejerce el muelle.
- ``L₀ \equiv \text{ Longitud natural del muelle}\,(m). (L₀ > 0)\longrightarrow \boldsymbol{L₀} = L₀ \,\boldsymbol{\hat{\imath}}``.
- ``-x \equiv \text{Elongación del muelle } (m)``. La elongación en esta primera fase puede ser negativa o nula; en el primer caso se debe preceder del signo negativo. Su valor absoluto es ``x`` ``\longrightarrow \boldsymbol{x} = -x\,\boldsymbol{\hat{\imath}}``.
- ``-x₀ \equiv \text{Elongación inicial del muelle}``. Es negativa porque el muelle se encuentra comprimido inicialmente. Su valor absoluto es ``x₀``. Puede tomar los valores: ``0 < x₀ \leq L₀``.
- ``\text{x₁} \equiv \text{Posición de } m₁``. Vale cero en toda esta fase ``(\mathbf{x₁} = 0\,\hat{\imath})``.
- ``\text{x₂} \equiv \text{Posición de } m₂ \, (\text{x₂} > 0)`` en toda esta fase ``(\text{\bf x₂} = \text{x₂}\,\boldsymbol{\hat{\imath}})``.
- ``\text{v₁} \equiv \text{Velocidad de } m₁``. Vale cero en toda la fase ``(\mathbf{v₁} = 0\,\boldsymbol{\hat{\imath}})``.
- ``\text{v₂} \equiv \text{Velocidad de } m₂`` en toda esta fase ``(\mathbf{v₂} = \pm \text{v₂}\,\boldsymbol{\hat{\imath}})``.
- ``\text{a₁} \equiv \text{Aceleración de } m₁``. Vale cero en toda la fase ``(\mathbf{a₁} = 0\,\boldsymbol{\hat{\imath}})``.
- ``\text{a₂} \equiv \text{Velocidad de } m₂`` en toda esta fase ``(\mathbf{a₂} = \pm \text{a₂}\,\boldsymbol{\hat{\imath}})``.
"""

# ╔═╡ aaaf469c-42e3-4d23-b57b-8793db22dbb0
md"#### Valores iniciales"

# ╔═╡ 7e519bb7-dc56-45df-9f44-e8d530b9a049
md"""
- ``\mathbf{x₁₀} = 0\,\boldsymbol{\hat{\imath}}\hspace{0.4em}(m) \longrightarrow \text{x₁₀} = 0\hspace{0.4em}(m)``.
- ``\mathbf{v₁₀} = 0\,\boldsymbol{\hat{\imath}}\hspace{0.4em}(m/s) \longrightarrow \text{v₁₀} = 0\hspace{0.4em}(m/s)``.
- ``\mathbf{x₂₀} = \boldsymbol{L₀} + \boldsymbol{x₀} \hspace{0.4em}(m) \longrightarrow \text{x₂₀} = L₀ - x₀ \hspace{0.4em}(m)``.
- ``\mathbf{v₂₀} = 0\,\boldsymbol{\hat{\imath}} \hspace{0.4em}(m/s) \longrightarrow v₂₀ = 0\hspace{0.4em}(m/s)``.\
"""

# ╔═╡ 21606d10-ebe1-419d-b463-7dd2cf20fe71
md"#### Aplicación de las leyes de Newton a la primera fase"

# ╔═╡ 930be4ad-2eee-47c7-90d5-e5ad2ccb4608
md"##### Sobre ``m₁``"

# ╔═╡ e482e2e5-7e90-4648-9d23-3c2b97ee8adb
md"""
- Sobre ``m₁`` actúan cuatro fuerzas. En el eje ``\text{x}`` actúa la fuerza elástica del muelle sobre la masa ``\boldsymbol{F}ₑ₁ = -kx\,\boldsymbol{\hat{\imath}}``, y en sentido opuesto, la reacción ``\boldsymbol{N}ₚ`` de la pared sobre ella. En el eje ``\text{y}`` actúan su peso ``m₁\boldsymbol{g}`` y la reacción del suelo ``\boldsymbol{N}ₛ₁``. Cuando la reacción de la pared se anule, esta masa estará libre para moverse; en este punto terminará esta primera fase.
"""

# ╔═╡ 83d12ae9-3c71-4210-80f9-39e592d7dfe5
md"$
\begin{gather}
\text{Eje x:}\hspace{1em} \boldsymbol{N}ₚ + \boldsymbol{F}ₑ₁ = 0 \hspace{1em}\longrightarrow\hspace{1em} Nₚ = kx\\
\hspace{1.45em}\text{Eje y:}\hspace{1em} \boldsymbol{N}ₛ₁ + m₁\boldsymbol{g}=0\hspace{1em}\longrightarrow\hspace{1em}Nₛ₁ = m₁ g
\end{gather}$
"

# ╔═╡ 83aa9ded-a4c4-47da-99db-053c67ee587b
md"##### Sobre ``m₂``"

# ╔═╡ f02e84e1-7ddb-4b46-b7ac-185ade9d9d83
md"""
- Sobre ``m₂`` actúan tres fuerzas. En el eje ``\text{x}`` actúa la fuerza elástica del muelle ``\boldsymbol{F}ₑ₂ = kx\,\boldsymbol{\hat{\imath}}``. En el eje ``\text{y}`` actúan su peso ``m₂\boldsymbol{g}``, y la reacción del suelo ``\,{N}ₛ₂``. En el eje ``\text{x}`` la fuerza elástica provoca una aceleración en ``m₂``.
"""

# ╔═╡ 6756649b-92e3-44e9-ab63-05db279b45df
md"$
\begin{gather}
\text{Eje x:}\hspace{1em} \boldsymbol{F}ₑ₂ = m₂\,\mathbf{a₂}\hspace{1em}\longrightarrow\hspace{1em}kx = m₂\,a₂ \\
\hspace{1.6em}\text{Eje y:}\hspace{1em} \boldsymbol{N}ₛ₂ + m₂\,\boldsymbol{g} =0\hspace{1em}\longrightarrow\hspace{1em}Nₛ₂ = m₂ \,g
\end{gather}$
"

# ╔═╡ 87744172-a35c-4360-a016-8f55ea9a67e3
md"#### Ecuación diferencial"

# ╔═╡ 50f66526-988c-41e1-93c4-b8899f90f3e0
md"""
``
\boldsymbol{Fₑ₂}(t) = m₂, \mathbf{a₂}(t)
``
"""

# ╔═╡ 07db6630-e83e-4f14-8fce-65fdb00c2153
md"""
``
\dfrac{d^2\mathbf{x₂}(t)}{dt^2} - \dfrac{k}{m₂} x(t) = 0
``
"""

# ╔═╡ bb104dc8-fd20-45fe-b0b3-044038da1cd3
md"Escribiendo la ecuación mediante escalares por tratarse de un movimiento unidimensional"

# ╔═╡ f446b922-4640-44f8-acba-388089d7331a
md"""
``
\dfrac{d^2 x₂(t)}{dt^2} + \dfrac{k}{m₂} (x₂(t) - L₀) = 0
``
"""

# ╔═╡ 18d770e4-0aca-4370-aceb-2f995aa64bc9
md"#### Cálculo analítico"

# ╔═╡ 0691719a-f48e-4c62-b1ac-47d906bdb5d1
md"""
Cambio de variable: ``X(t) = x₂(t) - L₀``
"""

# ╔═╡ ce31b020-9e41-4885-8085-48c479556735
md"""
``
\dfrac{d^2 X(t)}{dt^2} + \dfrac{k}{m₂} X(t) = 0
``
"""

# ╔═╡ f7d6c8f3-aa8e-4637-a2ee-7b0512465119
md"Solución de la ecuación diferencial anterior"

# ╔═╡ 4581a34a-225a-4960-bf84-efbf63bf0039
md"""
``
X = A \cos(ω₁ t + \delta);
\hspace{1em}
ω₁ = \sqrt{k/m₂}
``
"""

# ╔═╡ fd0bada9-6c77-4725-9835-99d64db36302
md"##### Solución general de la ecuación diferencial"

# ╔═╡ 54c7087b-e18b-41ff-9ed1-cfb231bf8dbc
md"$
\begin{gather}
\text{x₂}(t) &= L₀ + A \cos(ω₁ t + δ) \\
\text{v₂}(t) &= -A\,ω₁ \sin(ω₁ t + δ)
\end{gather}$
"

# ╔═╡ 14969615-02d6-48be-9d42-30ba0e2175cb
md"Condiciones iniciales:"

# ╔═╡ 76b1fe59-534d-4d45-86b0-f758f2a34ed9
md"""
* ``t = 0: \text{ x₂} = L₀ - x₀``
``\hspace{2.1em}A = \dfrac{-x₀}{\cos(\delta)}``;
``\hspace{0.5em}``
Para que la amplitud sea positiva:
``\hspace{0.5em}``
``\cos(δ) < 0``
"""

# ╔═╡ e937e066-2b5c-4585-99cc-987e72fd5f8f
md"""
* ``t = 0: \text{ v₂} = 0``
``\hspace{2.1em} \sin(δ) = 0 \hspace{1em}\longrightarrow\hspace{1em} δ = π``
"""

# ╔═╡ f21e13e0-6a37-4c69-8cf3-f3bdf6da0792
md"##### Solución particular de la ecuación diferencial"

# ╔═╡ 56710ca1-077e-449a-b086-dc1e69d7f558
md"$
\begin{gather}
\hspace{-7.2em}\text{x₁}(t) = 0\\
\hspace{-7.2em}\text{v₁}(t) = 0\\
\hspace{-7.2em}\text{a₁}(t) = 0\\
\hspace{0.8em}\text{x₂}(t) = L₀ + x₀ \cos(ω₁ t + π) \\
\hspace{-0.8em}\text{v₂}(t) = x₀\, ω₁ \sin(ω₁ t + π)\\
\hspace{-2.0em}\text{a₂}(t) = x₀\,ω₁^2 \cos(ωt)
\end{gather}$
``\hspace{13.7em}``con
``\hspace{1em}ω₁ = \sqrt{k/m₂}``
$"

# ╔═╡ d7725e84-5f41-4a15-8d0c-a79208d81023
md"La solución se puede simplificar"

# ╔═╡ 05a9576b-dc56-4621-a00a-9ab1033b8344
md"$
\begin{gather}
\hspace{1.6em}\text{x₂}(t) = L₀ - x₀ \cos(ω₁ t) \\
\text{v₂}(t) = x₀\, ω₁ \sin(ω₁ t)\\
\hspace{0.6em}\text{a₂}(t) = x₀\, ω₁^2 \cos(ω₁ t)
\end{gather}$
``\hspace{15em}``con
``\hspace{1em}ω₁ = \sqrt{k/m₂}``
$"

# ╔═╡ 2d7d3def-2cf1-4350-aeea-2c2b1db507cc
md"Esta fase termina cuando ``\text{x₂}(t) = L₀``"

# ╔═╡ 7a884cb3-5f42-4c34-a5b8-81ff7553e495
md"""
``L₀ - x₀ \cos(ω₁ t) = L₀
\hspace{1em}\longrightarrow\hspace{1em} \cos(ω₁ t) = 0
\hspace{1em}\longrightarrow\hspace{1em} t₁ = \dfrac{π}{2\,ω₁}``
"""

# ╔═╡ c1797c7a-1b0a-436f-b7fd-841d9b5e3bc5
md"##### Valores iniciales de m₁ y m₂"

# ╔═╡ 7e6e1f69-fff3-495d-b73e-43e551d9b9a6
md"$
\begin{gather}
\text{x₁₀} = 0\\
\text{v₁₀} = 0\\
\text{a₁₀} = 0\\
\hspace{2.7em}\text{x₂₀} = L₀ - x₀\\
\text{v₂₀} = 0\\
\hspace{1.6em}\text{a₁₀} = x₀\,\omega₁^2
\end{gather}$
"

# ╔═╡ 4fea284e-df13-4110-b587-25c2ab7c413e
md"##### Valores finales de m₁ y m₂"

# ╔═╡ 601c518b-b860-438a-8777-f370b1057488
md"$
\begin{gather}
\text{x₁ₚ} = 0\\
\text{v₁ₚ} = 0\\
\text{a₁ₚ} = 0\\
\hspace{0.5em}\text{x₂ₚ} = L₀\\
\hspace{1.0em}\text{v₂ₚ} = x₀\,ω₁\\
\text{a₂ₚ} = 0
\end{gather}$
"

# ╔═╡ af858578-9c8e-4b0f-ba80-6ddcc030e840
md"##### Centro de masas en función del tiempo"

# ╔═╡ 34b7d3fc-5467-4ec4-bc18-33f69f8779b0
md"$
\begin{gather}
\hspace{2.2em}\text{X}(t) = \dfrac{m₂}{m₁ + m₂} (L₀ - x₀\cos(ω₁ t)) = \dfrac{m₂}{m₁+m₂} \text{x₂}(t)\\
\text{V}(t) = \dfrac{m₂}{m₁ + m₂} x₀\,ω₁\sin(ω₁ t) = \dfrac{m₂}{m₁+m₂} \text{v₂}(t)\\
\hspace{0.5em}\text{A}(t) = \dfrac{m₂}{m₁ + m₂} x₀\,ω₁^2\cos(ω₁ t) = \dfrac{m₂}{m₁+m₂} \text{a₂}(t)
\end{gather}$
``\hspace{9.0em}``con
``\hspace{1em}ω₁ = \sqrt{k/m₂}``
$"

# ╔═╡ 56311590-1ce5-4316-92f4-6872cbcab85b
md"##### Valores iniciales del centro de masas"

# ╔═╡ 7da08f46-7865-47b1-86d0-07e74e32ed6f
md"$
\begin{gather}
\hspace{2.2em}\text{X₀} = \dfrac{m₂}{m₁ + m₂} (L₀ - x₀)\\
\text{V₀} = 0\\
\hspace{0.5em}\text{A₀} = \dfrac{m₂}{m₁ + m₂} x₀\,ω₁^2
\end{gather}$
``\hspace{15em}``con
``\hspace{1em}ω₁ = \sqrt{k/m₂}``
"

# ╔═╡ 554258ce-0ccb-47bc-a629-9a2811cfbf6f
md"##### Valores finales del centro de masas"

# ╔═╡ 29f24328-0805-4360-9222-d859987dcfb8
md"$
\begin{gather}
\hspace{-0.8em}\text{Xₚ} = \dfrac{m₂}{m₁ + m₂} L₀\\
\text{Vₚ} = \dfrac{m₂}{m₁ + m₂} x₀\,ω₁\\
\hspace{-4.9em}\text{Aₚ} = 0
\end{gather}$
``\hspace{15.2em}``con
``\hspace{1em}ω₁ = \sqrt{k/m₂}``
"

# ╔═╡ 891335e5-a97d-42a2-a9ee-cde6bd183f9c
md"##### Energía mecánica (fase 1)"

# ╔═╡ 606f1855-6040-4f5d-b5b1-3db3bed2c658
md"""
``
\hspace{16em}E = \dfrac{1}{2} k\,x₀^2 = \dfrac{1}{2} m₂\,ω₁^2;
\hspace{2em}\text{con } ω₁ = \sqrt{k/m₂}
``
"""

# ╔═╡ f74a28f8-cb6a-4fd9-8d6c-0cce6eb72611
md"### Segunda fase"

# ╔═╡ b93586da-4555-4b04-8685-bfdbeebb6da2
md"#### Parámetros"

# ╔═╡ bfee0f02-0c8e-463b-b5d7-f025a224acff
md"""
- ``m₁ \equiv \text{ Masa situada contra la pared (a la izquierda del sistema) } (kg)``. En esta segunda fase se separa de la pared.
- ``m₂ \equiv \text{ Masa situada a cierta distancia de la pared (a la derecha del sistema) } (kg)``. Esta masa se mueve sujeta a la fuerza que le ejerce el muelle.
- ``L₀ \equiv \text{ Longitud natural del muelle}\,(m). (L₀ > 0)\longrightarrow \boldsymbol{L₀} = L₀ \,\boldsymbol{\hat{\imath}}``.
- ``x \equiv \text{Elongación del muelle } (m)``. La elongación en esta segunda fase puede ser negativa, nula o positiva; en el primer caso se debe preceder del signo negativo. Su valor absoluto es ``x``. La elongación inicial del muelle en esta fase es cero ``\longrightarrow \boldsymbol{x} = \pm x\,\boldsymbol{\hat{\imath}}``.
- ``\text{x₁} \equiv \text{Posición de } m₁``. ``(\text{x₁} \geq 0 (\mathbf{x₁} = \pm\text{x₁}\,\hat{\imath})``.
- ``\text{x₂} \equiv \text{Posición de } m₂ \, (\text{x₂} > 0)`` en toda esta fase ``(\text{\bf x₂} = \text{x₂}\,\boldsymbol{\hat{\imath}})``.
- ``\text{v₁} \equiv \text{Velocidad de } m₁`` en toda esta fase ``(\mathbf{v₁} = \pm \text{v₁}\,\boldsymbol{\hat{\imath}})``.
- ``\text{v₂} \equiv \text{Velocidad de } m₂`` en toda esta fase ``(\mathbf{v₂} = \pm \text{v₂}\,\boldsymbol{\hat{\imath}})``.
- ``\text{a₁} \equiv \text{Aceleración de } m₁`` en toda esta fase ``(\mathbf{a₁} = \pm \text{a₁}\,\boldsymbol{\hat{\imath}})``.
- ``\text{a₂} \equiv \text{Aceleración de } m₂`` en toda esta fase ``(\mathbf{a₂} = \pm \text{a₂}\,\boldsymbol{\hat{\imath}})``.
"""

# ╔═╡ 2ae54de3-9506-49ac-b238-5d944089b7b7
md"#### Valores iniciales"

# ╔═╡ 5f523482-86a2-4f36-9f8e-ef3a4d72e53a
md"""
- ``\mathbf{x₁₀} = 0\,\boldsymbol{\hat{\imath}}\hspace{0.4em}(m) \longrightarrow \text{x₁₀} = 0\hspace{0.4em}(m)``.
- ``\mathbf{v₁₀} = 0\,\boldsymbol{\hat{\imath}}\hspace{0.4em}(m/s) \longrightarrow \text{v₁₀} = 0\hspace{0.4em}(m/s)``.
- ``\mathbf{x₂₀} = \boldsymbol{L₀} \hspace{0.4em}(m) \longrightarrow \text{x₂₀} = L₀ \hspace{0.4em}(m)``.
- ``\mathbf{v₂₀} = x₀\,ω₁\,\boldsymbol{\hat{\imath}} \hspace{0.4em}(m/s) \longrightarrow v₂₀ = x₀\,ω₁\hspace{0.4em}(m/s)``.\
"""

# ╔═╡ ebe29b79-59c2-4cfa-b964-9d978345bd91
md"#### Aplicación de las leyes de Newton a la segunda fase"

# ╔═╡ ea584508-2b28-4942-80a4-1fccd9f4e954
md"##### Sobre ``m₁``"

# ╔═╡ a22f3e68-eed3-4182-b169-046464a87172
md"""
- Sobre ``m₁`` actúan ahora tres fuerzas. En el eje ``\text{x}`` actúa la fuerza elástica del muelle sobre la masa ``\boldsymbol{F}ₑ₁ = -kx\,\boldsymbol{\hat{\imath}}``, pero ya no hay reacción ``\boldsymbol{N}ₚ`` de la pared sobre ella. En el eje ``\text{y}`` actúan su peso ``m₁\boldsymbol{g}`` y la reacción del suelo ``\boldsymbol{N}ₛ₁``. Cuando la reacción de la pared se anule, esta masa estará libre para moverse; en este punto terminará esta primera fase.
"""

# ╔═╡ 7874884f-2cb1-4c71-a757-4ae5cff6b5c0
md"$
\begin{gather}
\text{Eje x:}\hspace{1em} \boldsymbol{F}ₑ₁ = m₁ \mathbf{a₁} \hspace{1em}\longrightarrow\hspace{1em} kx = m₁ \text{a₁}\\
\hspace{1.48em}\text{Eje y:}\hspace{1em} \boldsymbol{N}ₛ₁ + m₁\boldsymbol{g}=0\hspace{1em}\longrightarrow\hspace{1em}Nₛ₁ = m₁ g
\end{gather}$
"

# ╔═╡ e912655a-f2bd-4d7f-aaea-2612c25a5b3e
md"##### Sobre ``m₂``"

# ╔═╡ b6319630-97cf-43d2-bb73-aedeeaa903a1
md"""
- Sobre ``m₂`` actúan tres fuerzas. En el eje ``\text{x}`` actúa la fuerza elástica del muelle ``\boldsymbol{F}ₑ₂ = kx\,\boldsymbol{\hat{\imath}}``. En el eje ``\text{y}`` actúan su peso ``m₂\boldsymbol{g}``, y la reacción del suelo ``\,{N}ₛ₂``. En el eje ``\text{x}`` la fuerza elástica provoca una aceleración en ``m₂``.
"""

# ╔═╡ dc48a744-10d9-4885-982f-ebca319f4aee
md"$
\begin{gather}
\text{Eje x:}\hspace{1em} \boldsymbol{F}ₑ₂ = m₂\,\mathbf{a₂}\hspace{1em}\longrightarrow\hspace{1em}-kx = m₂ \,\text{a₂} \\
\hspace{1.8em}\text{Eje y:}\hspace{1em} \boldsymbol{N}ₛ₂ + m₂\,\boldsymbol{g} =0\hspace{1em}\longrightarrow\hspace{1em}Nₛ₂ = m₂\, g
\end{gather}$
"

# ╔═╡ 3668c039-e55e-4bf3-959c-1613f9b8b481
md"#### Ecuaciones diferenciales acopladas"

# ╔═╡ afe05bcd-357d-4869-83fd-89a398504539
md"$
\begin{cases}
m₁\,\text{a₁}(t) = k\,x(t)\\
m₂\,\text{a₂}(t) = -k\,x(t)
\end{cases}$
``\hspace{14.8em}\text{Vemos que: } m₁\,\text{a₁} = -m₂\,\text{a₂}``
"

# ╔═╡ 2a30b827-a4ac-4f15-914f-b2941c6ac10a
md"$
\begin{cases}
\dfrac{d^2\text{x₁}(t)}{dt^2} = \dfrac{k}{m₁}\,x(t)\\
\dfrac{d^2\text{x₂}(t)}{dt^2} = -\dfrac{k}{m₂}\,x(t)
\end{cases}$"

# ╔═╡ ec6511a6-baa1-4959-8f45-4f1deb9d49fd
md"###### Ecuaciones diferenciales acopladas que resolveremos mediante `DifferentialEquations.jl` en `Julia`"

# ╔═╡ 958c56e3-26ef-4cf2-b3e8-ce7ff5b9d7e6
md"$
\begin{cases}
\dfrac{d^2\text{x₁}(t)}{dt^2} = \dfrac{k}{m₁}\,(\text{x₂} - \text{x₁} - L₀)\\
\dfrac{d^2\text{x₂}(t)}{dt^2} = -\dfrac{k}{m₂}\,(\text{x₂} - \text{x₁} - L₀)
\end{cases}$"

# ╔═╡ ac175b42-b5a9-4462-805b-1de215710578
md"#### Cálculo analítico"

# ╔═╡ 6ec06427-b723-4fc0-8bec-8a92b74417c0
md"Para analizar analíticamente este sistema de ecuaciones, restamos la primera de la segunda"

# ╔═╡ 0a3d31d0-5c25-4d18-b176-22f1eee61cd5
md"""
```math
\dfrac{d^2\text{x₂}}{dt^2} - \dfrac{d^2\text{x₁}}{dt^2}
= \dfrac{k}{m₂}(\text{x₂}-\text{x₁}-L₀)-\dfrac{k}{m₁}(\text{x₂}-\text{x₁}-L₀)
```
"""

# ╔═╡ 95b48e57-9b22-442d-9f79-24e3ab148bf3
md"""
```math
\dfrac{d^2\text{x₂₁}}{dt^2}
= -k\left(\dfrac{1}{m₂}+\dfrac{1}{m₁}\right)\,(\text{x₂₁}-L₀)
```
"""

# ╔═╡ f8e14404-bdc8-4f48-b278-fbcf6a4566dd
md"""
```math
\dfrac{d^2\text{x₂₁}}{dt^2}
= -\dfrac{k}{\mu}\,(\text{x₂₁}-L₀)
```
"""

# ╔═╡ 5bcb6eec-59b9-4096-a95b-c136281baad1
md"donde ``μ`` es la masa reducida del sistema"

# ╔═╡ 3a788f13-59b8-48a3-ab9f-48b96d2a3486
md"""
Cambio de variable: ``\text{X₂₁}(t) = \text{x₂₁}(t) - L₀``
"""

# ╔═╡ be85ce71-110a-4dd5-bbfd-52f19784d2e0
md"""
```math
\dfrac{d^2 X₂₁(t)}{dt^2} + \dfrac{k}{μ} X₂₁(t) = 0
```
"""

# ╔═╡ 381e74e7-4754-4ad7-8625-520c3674c666
md"Solución de la ecuación diferencial anterior"

# ╔═╡ 1ba5ad0c-be87-4629-9c19-f81d54cad92b
md"""
```math
\text{X₂₁} = A \cos(ω₂ t + \delta);
\hspace{1em}
ω₂ = \sqrt{k/μ}
```
"""

# ╔═╡ 064cd216-775e-4995-b13c-ed00c8e56b51
md"""
```math
\text{x₂} - \text{x₁} = L₀ + A \cos(ω₂\,t + \delta)
```
"""

# ╔═╡ c296ceca-501c-4cdb-bc76-6a9f796d38cb
md"##### Solución general de la ecuación diferencial"

# ╔═╡ e7424c03-0db4-447e-8f17-687a6f01347c
md"$
\begin{gather}
\text{x₂}(t) - \text{x₁}(t) = L₀ + A \cos(ω₂\,t + δ) \\
\hspace{-0.5em}\text{v₂}(t) - \text{v₁}(t) = -A\,ω₂ \sin(ω₂\,t + δ)
\end{gather}$
"

# ╔═╡ c159a02c-2384-4c11-a1e5-03438de7cda6
md"Condiciones iniciales:"

# ╔═╡ 75494f1e-a3f7-4790-b364-f58fa8910e42
md"""
* ``t = \dfrac{π}{2\,ω₁}: \text{x₂₀} - \text{x₁₀} = L₀ - 0 = L₀``
``\hspace{2.1em}L₀ = L₀ + A \cos(δ)``
``\hspace{1em}\longrightarrow\hspace{1em}``
``\cos(\delta) = 0 \hspace{1em}\longrightarrow\hspace{1em} δ = \dfrac{π}{2}, \dfrac{3π}{2}, \cdots``
"""

# ╔═╡ 299f987d-fa02-4fe7-901d-d2c40d8e785a
md"""
* ``t = \dfrac{π}{2\,ω₁}: \text{v₂₀} - \text{v₁₀} = x₀\,ω₁ - 0 = x₀\,ω₁``
``\hspace{2.1em}x₀\,ω₁ = -A\,ω₂ \sin(\delta)``
``\hspace{1em}\longrightarrow\hspace{1em}``
``A = x₀ \dfrac{ω₁}{ω₂\,\sin(δ)} > 0\hspace{1em}\longrightarrow\hspace{1em} δ = \dfrac{\pi}{2}``
"""

# ╔═╡ 88630861-b840-44cf-8230-39436b2888fc
md"##### Solución particular de la ecuación diferencial"

# ╔═╡ 2736abc6-c8fa-42f9-a7f1-8c4797dea825
md"$
\begin{gather}
\hspace{-5.6em}\text{x₂}(t) - \text{x₁}(t) = L₀ + x₀ \dfrac{ω₁}{ω₂} \cos\left(ω₂\,t + \dfrac{π}{2}\right)\\
\hspace{-7.2em}\text{v₂}(t) - \text{v₁}(t) = -x₀\,ω₁ \sin\left(ω₂\,t + \dfrac{π}{2}\right)\\
\hspace{-6.0em}\text{a₂}(t) - \text{a₁}(t) = -x₀\,ω₁\,ω₂ \cos\left(ω₂\,t + \dfrac{π}{2}\right)\\
\end{gather}$
``\hspace{7.4em} \text{ con } \hspace{1em}\dfrac{ω₁}{ω₂} = \sqrt{\dfrac{m₁}{m₁+m₂}}
\hspace{2em}\text{y}\hspace{2em} ω₁\,ω₂ = ω₁^2 \,\sqrt{\dfrac{m₁+m₂}{m₁}}``
"

# ╔═╡ 2dc09545-4c46-4833-8d76-45e6383fad7d
md"La solución se puede simplificar"

# ╔═╡ 173325de-9ac5-48a1-912c-3806212a1796
md"$
\begin{gather}
\hspace{-5.6em}\text{x₂}(t) - \text{x₁}(t) = L₀-x₀\,\dfrac{ω₁}{ω₂}\,\sin(ω₂\,t)\\
\hspace{-7.2em}\text{v₂}(t) - \text{v₁}(t) = -x₀\,ω₁\,\cos(ω₂\,t)\\
\hspace{-7.0em}\text{a₂}(t) - \text{a₁}(t) = x₀\,ω₁\,ω₂\,\sin(ω₂\,t)\\
\end{gather}$
``\hspace{7.4em} \text{ con } \hspace{1em}\dfrac{ω₁}{ω₂} = \sqrt{\dfrac{m₁}{m₁+m₂}}
\hspace{2em}\text{y}\hspace{2em} ω₁\,ω₂ = \dfrac{k}{m₂} \,\sqrt{\dfrac{m₁+m₂}{m₁}}``
"

# ╔═╡ 052bf6ad-9fcc-4037-9115-dcbb2f19d99b
md"##### Aceleración del centro de masas"

# ╔═╡ 84b4fc0b-ee8d-43d0-acb9-e1b3d8146224
md"""
```math
\text{A}(t)= \dfrac{m₁\,\text{a₁} + m₂\,\text{a₂}}{m₁ + m₂} = \dfrac{m₁\,\dfrac{kx}{m₁} + m₂\,\dfrac{-kx}{m₂}}{m₁ + m₂} = 0
```
"""

# ╔═╡ 17f19d51-a5a3-497f-ba26-68cdc7ff9834
md"##### Velocidad del centro de masas"

# ╔═╡ 8fbdd5bb-ba8c-4ba9-a90e-ae804c50d61b
md"El centro de masas no tiene aceleración. Tendrá velocidad constante. Por tanto, la velocidad del centro de masas en esta fase será la última de la fase anterior"

# ╔═╡ fb43062f-fed6-4120-bc4b-2af832c16387
md"""
```math
V(t) = \dfrac{m₂}{m₁+m₂} x₀\,ω₁
```
"""

# ╔═╡ 40b922a7-5155-4cf6-a50c-4597a32c2bf4
md"##### Posición del centro de masas"

# ╔═╡ 26dffbe9-d8c4-402b-9913-7172c0458b8e
md"""
```math
\text{X}(t) = X₀ + \text{V}(t)\,t
```
"""

# ╔═╡ 0c4207fd-3b47-4042-b821-ca24917145ed
md"""
```math
\text{X}(t) = \dfrac{m₂}{m₁+m₂} L₀ + \dfrac{m₂}{m₁+m₂} x₀\,ω₁\,t
```
"""

# ╔═╡ 9c9f9325-e5e1-4ebb-9402-fdbf49de8e32
md"""
```math
\text{X}(t) = \dfrac{m₂}{m₁+m₂} (L₀ + x₀\,ω₁\,t)
```
"""

# ╔═╡ 7f4dee6a-7e72-4531-9c42-d02fd8d3bf12
md"### Código fuente general"

# ╔═╡ 924f3ca7-a169-4297-846c-e88607d6fd54
md"#### Carga de paquetes"

# ╔═╡ 28c31f91-e55f-4d3b-b27d-a0ea722a7447
md"#### Botón reset"

# ╔═╡ 28a797d4-7e36-4283-9daa-aefae713c4a4
@bind reset Button("Reset")

# ╔═╡ 5a5bcb8d-f0fb-406f-acd4-680edc09c1df
md"##### Parámetros fijos"

# ╔═╡ 3dccada0-19c5-4dad-9616-f73a9d00d86d
md"###### L₀: Longitud natural del muelle en $m$"

# ╔═╡ e04a6c73-e3cf-48ab-908f-fd16798062a0
md"###### x₀: Elongación inicial del muelle en valor absoluto"

# ╔═╡ 5486689c-ab65-439f-98ce-e358e6ba24f2
md"###### k: Constante elástica del muelle en $N/m$"

# ╔═╡ 0c8af1d1-381d-4166-925b-246ba9e72d08
md"###### m₁: Masa más cercana a la pared"

# ╔═╡ f12c8235-b240-494d-91ee-0b264ba9a4c0
md"###### m₂: Masa más alejada de la pared"

# ╔═╡ 7807081f-c329-4e26-b4d2-0d2848d2bc0e
md"###### μ: Masa reducida del sistema"

# ╔═╡ 7d3ba095-977d-4108-857f-039dc50b022c
md"###### ω₁: Frecuencia angular en la primera fase"

# ╔═╡ b0534f3d-6481-4e3d-aa48-0888eeceab23
md"###### ω₂: Frecuencia angular en la segunda fase"

# ╔═╡ b0142e00-5f55-4e5f-bd49-962b862f6312
md"###### t₁₀: Tiempo inicial de la primera fase"

# ╔═╡ 236e9c09-571e-47a9-92c1-fa9ac4e6965e
md"###### t₁ₚ: Tiempo final de la primera fase"

# ╔═╡ 21eaf4af-2f31-4f2a-8cbf-c3719c1db356
md"###### t₂₀: Tiempo inicial de la segunda fase"

# ╔═╡ d952b13c-ef40-4f46-bfcd-883f951faf74
md"###### t₂ₚ: Tiempo final de la segunda fase"

# ╔═╡ c3ab23be-07f4-4429-834a-a612ae221a4a
begin
	L₀ = 0.25
	x₀ = 0.10
	k = 1.0
	m₁ = 1.0
	m₂ = 1.0
	μ = 1/m₁ + 1/m₂
	ω₁ = sqrt(k/m₁)
	ω₂ = sqrt(k/μ)
	t₀₁ = 0.0
	tₚ₁ = π/(2*ω₁)
	t₀₂ = tₚ₁
	tₚ₂ = tₚ₁ + 1 * 2*π/ω₂
	x₁₀₁ = 0.0
	x₂₀₁ = L₀ - x₀
	v₁₀₁ = 0.0
	v₂₀₁ = 0.0
	x₁ₚ₁ = x₁₀₁
	x₂ₚ₁ = L₀
	v₁ₚ₁ = v₁₀₁
	v₂ₚ₁ = x₀ * ω₁
	x₁₀₂ = x₁ₚ₁
	x₂₀₂ = x₂ₚ₁
	v₁₀₂ = v₁ₚ₁
	v₂₀₂ = v₂ₚ₁
end

# ╔═╡ c677f634-3133-44ce-bc5e-04c2d516f742
md"## Código fuente (fase 1)"

# ╔═╡ 5e3d0230-26be-42bb-a057-4784be1263f5
md"### Funciones auxiliares"

# ╔═╡ 065d6b1c-67de-4070-9c1e-774553b49168
md"#### Extracción de resultados en forma de arrays"

# ╔═╡ e0e48f7a-dcd7-4fb6-ae6a-1f3635618cfc
function valores_ODE_fase1(sol; dt=0.001)
	idx = sol.t[1]:dt:sol.t[end]
	return (collect(idx),
	        map(value -> value[1], sol.(idx)),
	        map(value -> value[2], sol.(idx)))
end

# ╔═╡ e660648a-8a93-4efc-9e58-5207307263be
md"### Planteamiento de la ecuación diferencial"

# ╔═╡ 18fc5377-34ae-4782-bf7e-62b741286044
md"#### Definición de la función ODE en `Julia`"

# ╔═╡ 57432286-84f8-4c66-84ba-2c03c5769d5e
function fase1(du, u, p, t)
	dx₂ = du
	x₂ = u
	L₀, x₀, k, m₂ = p
	-(k/m₂)*(x₂ - L₀)
end

# ╔═╡ 83deae0a-a7ce-400c-8dbe-07bce3d8e8f8
md"#### Parámetros que se le pasan a la función ODE"

# ╔═╡ ebb2888e-4561-44e3-998d-60241bca7c2d
p1 = [L₀, x₀, k, m₂]

# ╔═╡ 8e45f531-a6ed-48c7-a21c-c218dae6a9f3
md"#### Asignación de variables"

# ╔═╡ 433dd48b-cd87-48ef-86d6-8e32959c2fbf
md"##### ``du = dx₂ = v₂\,\equiv`` Velocidad de la masa más alejada de la pared ($m₂$) en $m/s$"

# ╔═╡ a31afff7-a77f-4b36-b576-7ef3f9f33380
md"##### ``u = x₂\,\equiv`` Posición de la masa más alejada de la pared ($m₂$) en $m$"

# ╔═╡ b97fb2e2-9f61-4555-b50f-94c5c69c5b71
md"#### Tiempo en $s$"

# ╔═╡ c2f0923e-58b9-4b2f-8c8d-f6174e6e4523
tspan1 = (t₀₁, tₚ₁)

# ╔═╡ 0fb3aebe-affb-4849-a987-5fe11c900fd9
md"#### Definición de la ODE de segundo orden"

# ╔═╡ d6eacbd0-7e55-4c8d-829f-54a4f47051b9
fase1_prob = SecondOrderODEProblem(fase1, v₂₀₁, x₂₀₁, tspan1, p1)

# ╔═╡ f475a8f0-d55b-48aa-b5e9-fcfe358add4e
md"##### Solución del problema"

# ╔═╡ 8d131541-edc6-4cd4-b35a-4bb4981dbee6
fase1_sol = solve(fase1_prob, abstol=1e-12, reltol=1e-12)

# ╔═╡ 03b0c959-9d0d-46e4-96c2-972a0771922f
md"###### Extracción de valores del resultado en forma de arrays"

# ╔═╡ 6a887098-abd4-41ba-9e50-de4b7eae2fc9
ODE_fase1_t, ODE_fase1_v₂, ODE_fase1_x₂ = valores_ODE_fase1(fase1_sol, dt=0.001)

# ╔═╡ 369d58f4-6c8d-48af-b172-751e89f4ceea
md"###### Tamaño de los vectores de datos en la fase 1"

# ╔═╡ d7ab580c-061f-4369-babf-1c7608d24073
fase1_length = length(ODE_fase1_t)

# ╔═╡ f76e194d-65da-48bc-ac2c-c213473f79e6
md"###### Posición ``\text{x₁}(t)``"

# ╔═╡ 0b53878d-9a5d-4b96-a2eb-5a2bdc1e4027
ODE_fase1_x₁ = zeros(fase1_length)

# ╔═╡ 1cfc9bbd-05cd-4d31-90e9-f8ba2858fb23
md"###### Posición ``\text{x₂}\,(t)``, calculada analíticamente"

# ╔═╡ 073cd0dc-4932-455b-9f88-a350e2e9695a
fase1_analitica_x₂ = L₀ .- x₀ .* cos.(ω₁ .* ODE_fase1_t)

# ╔═╡ 8b369410-4589-4878-996d-9cf034e1281f
md"###### Velocidad ``\text{v₁}\,(t)``"

# ╔═╡ 1552fbf4-a1c3-4daa-9ead-387c11dfa23a
ODE_fase1_v₁ = ODE_fase1_x₁

# ╔═╡ 372f8d9a-7cb1-4ff6-8cbe-b81c57872feb
md"###### Velocidad ``\text{v₂}(t)``, calculada analíticamente"

# ╔═╡ 7b531ed8-44b4-4a60-b182-24b1a0f5dc01
fase1_analitica_v₂ = x₀ .* ω₁ .* sin.(ω₁ .* ODE_fase1_t)

# ╔═╡ 80b47e79-e777-4d6c-9ef7-ffdea8fcf1ee
md"###### Aceleración ``\text{a₁}\,(t)``"

# ╔═╡ 0d733fd6-27cd-4a54-bced-9f579a7c3ef3
ODE_fase1_a₁ = ODE_fase1_x₁

# ╔═╡ 9175746e-286f-4856-ad16-f64d479f7842
md"###### Aceleración ``\text{a₂}\,(t)``, calculada analíticamente"

# ╔═╡ ff64e888-2731-409c-86f9-2d10418375a7
fase1_analitica_a₂ = x₀ .* ω₁^2 .* cos.(ω₁ .* ODE_fase1_t)

# ╔═╡ db21bff9-51ad-490e-ae9f-cc2deea8dfeb
md"###### Aceleración ``\text{a₂}\,(t)``"

# ╔═╡ 61cf579c-3807-4234-8a40-99ccd1cd54f1
ODE_fase1_a₂ = (k/m₂) .* (L₀ .- ODE_fase1_x₂)

# ╔═╡ a89390ae-cbdd-4144-8f0c-ed30e8c8f8d7
md"###### Posición del centro de masas del sistema ``\text{Xcm}\,(t)``"

# ╔═╡ f3b6b1a2-b267-4035-9267-b8830318de17
ODE_fase1_Xcm = (m₂/(m₁ + m₂)) .* ODE_fase1_x₂

# ╔═╡ 68cd7717-d64f-4298-9ff1-1f4ed61b5b86
md"###### Velocidad del centro de masas del sistema ``\text{Vcm}\,(t)``"

# ╔═╡ d1054a56-c567-4f82-aaac-5218e09e71cd
ODE_fase1_Vcm = (m₂/(m₁ + m₂)) .* ODE_fase1_v₂

# ╔═╡ 5dd36860-bd12-41e3-89a2-0fec5b2ad710
md"###### Aceleración del centro de masas del sistema ``\text{Acm}\,(t)``"

# ╔═╡ 729d72ab-51f7-410f-852a-3d85cd39ff4c
ODE_fase1_Acm = (m₂/(m₁ + m₂)) .* ODE_fase1_a₂

# ╔═╡ e556f528-b5b7-437b-b325-59af944e6697
md"###### Elongación del muelle ``x(t)``"

# ╔═╡ 89cf5e1f-d657-45ac-82ac-edb72c1aeb3b
ODE_fase1_elong = L₀ .- ODE_fase1_x₂

# ╔═╡ b35f2726-e632-466a-8f5e-e514c4c53697
md"###### Energía cinética del sistema ``Ec(t)``"

# ╔═╡ 92899b96-f19e-446f-96fa-1fc2c2d6e04d
ODE_fase1_Ec = 0.5 .* (m₁ .* ODE_fase1_v₁.^2 + m₂ .* ODE_fase1_v₂.^2)

# ╔═╡ 0cd13e27-cffa-43c0-b5de-cbc0158845e8
md"###### Energía potencial del muelle ``Ep(t)``"

# ╔═╡ 244fa152-13ab-45f5-b7f4-a66067fd53d0
ODE_fase1_Ep = 0.5 .* k .* ODE_fase1_elong.^2

# ╔═╡ ac0ca9aa-32d9-4e71-a3cf-a1f95e253a51
md"###### Energía total ``E(t)``, formado por masas y muelle"

# ╔═╡ 13c7f234-0abf-4886-8989-48e8f2ce4b5c
ODE_fase1_E = ODE_fase1_Ec .+ ODE_fase1_Ep

# ╔═╡ cb0ffc15-5abd-4bfa-ae0f-c05769f2f386
md"---"

# ╔═╡ 2a42e1e2-5123-4888-8281-164d44229425
md"#### Gráficas (fase 1)"

# ╔═╡ b9b8faee-e7f6-4bcf-b4a7-e6cefe5e2341
begin
	plot_fase1_x = plot(title = "Posición de m₂, del CM y elongación (fase 1)",
						 xaxis = "tiempo (s)", yaxis = "posición (m)",
						 legend_columns = 2, legend = :outerbottom);
	plot!(plot_fase1_x, ODE_fase1_t, ODE_fase1_x₁,
		  linecolor = :blue, alpha = 0.6, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "x₁(t)"
		 );
	plot!(plot_fase1_x, ODE_fase1_t, ODE_fase1_x₂,
		  linecolor = :red, alpha = 0.3, yguidefontcolor = :black,
		  linewidth = 8.0,
		  label = "x₂(t) ODE"
		 );
	plot!(plot_fase1_x,
		  ODE_fase1_t, fase1_analitica_x₂,
		  linecolor = :black,
		  linewidth = 2,
		  label = "x₂(t) = L₀ - x₀ω₁cos(ω₁t)"
		 );
	plot!(plot_fase1_x,
		  ODE_fase1_t, ODE_fase1_Xcm,
		  linecolor = :grey,
		  linewidth = 2,
		  label = "Xcm(t)"
		 )
end

# ╔═╡ 44d5951c-9562-4567-9dea-bd6425ba4676


# ╔═╡ 60db7902-2381-43e1-9a16-8d8ff33ce556
begin
	plot_fase1_v = plot(title = "Velocidad de m₂ y del CM (fase 1)",
					    xaxis = "tiempo (s)", yaxis = "velocidad (m/s)"
					   );
	plot!(plot_fase1_v, ODE_fase1_t, ODE_fase1_v₁,
		  linecolor = :blue, alpha = 0.6, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "v₁(t)"
		 );	
	plot!(plot_fase1_v, ODE_fase1_t, ODE_fase1_v₂,
		  linecolor = :red, alpha = 0.3, yguidefontcolor = :black,
		  linewidth = 8.0,
		  xaxis = "tiempo (s)", yaxis = "velocidad (m/s)",
		  label = "v₂(t) ODE"
		 );
	plot!(plot_fase1_v,
		  ODE_fase1_t, fase1_analitica_v₂,
		  linecolor = :black,
		  linewidth = 2,
		  label = "v₂(t) = x₀ ω₁ sin(ω₁ t)"
		 );
	plot!(plot_fase1_v,
		  ODE_fase1_t, ODE_fase1_Vcm,
		  linecolor = :gray,
		  linewidth = 2,	  
		  label = "Vcm(t)"
		 )
end

# ╔═╡ 0f89d7bf-75e9-485a-bd48-27358295c0b0
begin
	plot_fase1_a = plot(title = "Aceleración de m₂ y del CM (fase 1)",
					    xaxis = "tiempo (s)", yaxis = "aceleración (m/s²)");
	plot!(plot_fase1_a, ODE_fase1_t, ODE_fase1_a₁,
		  linecolor = :blue, alpha = 0.6, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "a₁(t)"
		 );
	plot!(plot_fase1_a, ODE_fase1_t, ODE_fase1_a₂,
		  linecolor = :red, alpha = 0.3, yguidefontcolor = :black,
		  linewidth = 8.0,
		  label = "a₂(t) ODE"
		 );
	plot!(plot_fase1_a,
		  ODE_fase1_t, fase1_analitica_a₂,
		  linecolor = :black,
		  linewidth = 2.0,
		  label = "a₂(t) = x₀ ω₁^2 sin(ω₁ t)"
		 )
	plot!(plot_fase1_a,
		  ODE_fase1_t, ODE_fase1_Acm,
		  linecolor = :gray,
		  linewidth = 2,
		  label = "Acm(t)"
		 )
end

# ╔═╡ f43f4f52-d53e-4b46-98a1-fb80cc4caa83
md"## Código fuente (fase 2)"

# ╔═╡ 2ffbe9cc-90fe-4c80-a1ca-94accea6f63a
md"### Funciones auxiliares"

# ╔═╡ b622ca7d-2121-4617-bc2c-05c51b465a65
md"#### Extracción de resultados en forma de arrays"

# ╔═╡ 171bf0fd-a545-4161-ae77-f3a23aa19233
function valores_ODE_fase2(sol; vars=(1,2,3,4), dt=0.001)
	idx = sol.t[1]:dt:sol.t[end]
	return (collect(idx),
	        map(value -> value[1], sol.(idx)),
	        map(value -> value[2], sol.(idx)),
	        map(value -> value[3], sol.(idx)),
	        map(value -> value[4], sol.(idx))
		   )
end

# ╔═╡ f021828f-e7d6-46a1-90cf-017345536bdd
md"### Planteamiento del sistema de EDO acopladas"

# ╔═╡ 5b712a32-13e7-4504-bb8a-90e92d6d246a
md"#### Definición de la función ODE en `Julia`"

# ╔═╡ 703907bb-a878-4823-87d1-ed00fb0ff0ed
function fase2(du, u, p, t)
	x₁, x₂, v₁, v₂ = u
	L₀, x₀, k, m₁, m₂, μ, ω₁, ω₂ = p
	du[1] = v₁
	du[2] = v₂
	du[3] = (k/m₁) * (x₂ - x₁ - L₀)
	du[4] = -(k/m₂) * (x₂ - x₁ - L₀)
end

# ╔═╡ 5c60e568-d358-42b6-aaa6-97fd1f51302a
md"#### Parámetros que se le pasan a la función ODE"

# ╔═╡ fbe27e04-9ff7-4b63-b91d-08d7513cce32
p2 = [L₀, x₀, k, m₁, m₂, μ, ω₁, ω₂]

# ╔═╡ 84b89a14-4b20-4525-8887-c4ad35e6d74a
md"#### Asignación de variables"

# ╔═╡ feb81d34-8870-4fd0-ae46-0704d7917a9d
md"##### ``dx₁ = v₁ \,\equiv`` Velocidad de la masa más cercana a la pared ($m₁$) en $m/s$"

# ╔═╡ 2859858a-9918-450b-8851-a5995645ba84
md"##### ``dx₂ = v₂\,\equiv`` Velocidad de la masa más alejada de la pared ($m₂$) en $m/s$"

# ╔═╡ 82f1437d-704a-4566-b70e-1061cbbf8b7a
md"##### ``x₁\,\equiv`` Posición de la masa más cercana a la pared ($m₁$) en $m$"

# ╔═╡ 3ebf307c-089e-4090-b875-d599902c63f7
md"##### ``x₂\,\equiv`` Posición de la masa más alejada de la pared ($m₂$) en $m$"

# ╔═╡ 2793e637-60f4-46c0-8de5-a699cf4554d0
md"#### Valores iniciales de posición y velocidad"

# ╔═╡ a73a864c-2603-4952-818c-3f48484b40cf
u_begin_fase2 = [x₁₀₂, x₂₀₂, v₁₀₂, v₂₀₂]

# ╔═╡ 0d68a866-5182-489b-bb0e-bf71e63bdb7e
md"#### Tiempo en $s$"

# ╔═╡ 169b2496-f2f6-4966-b971-b6615ea8e241
tspan2 = (t₀₂, tₚ₂)

# ╔═╡ b6f47d8f-c1c9-4a89-a665-1157daf728f7
md"#### Definición del sistema de ecuaciones diferenciales de segundo orden"

# ╔═╡ 6377237c-44cf-41f7-b8b8-b901ab4a1dd6
fase2_prob = ODEProblem(fase2, u_begin_fase2, tspan2, p2)

# ╔═╡ a2448f50-744e-443c-8ee1-d2730f7b9161
md"##### Solución del problema"

# ╔═╡ 849c2bce-edb3-4796-9414-c54625aa0977
fase2_sol = solve(fase2_prob, abstol=1e-12, reltol=1e-12)

# ╔═╡ e107ca77-2312-4fcd-a22a-58d524ac8aec
md"###### Extracción de valores del resultado en forma de arrays"

# ╔═╡ 8ca72f78-588d-4c69-a8e8-4fa16fb4b93f
ODE_fase2_t, ODE_fase2_x₁, ODE_fase2_x₂, ODE_fase2_v₁, ODE_fase2_v₂ = valores_ODE_fase2(fase2_sol, dt=0.001)

# ╔═╡ ab04a12a-b375-4038-ad6c-879033ffb3b6
ODE_fase1_x₂

# ╔═╡ 4fa8e88a-f332-4446-88f2-90abc693723f
md"###### Posición relativa ``\text{x₂}\,(t) - \text{x₁}\,(t)``"

# ╔═╡ 60f893c6-429a-4bbd-a534-b0d988afb5e4
ODE_fase2_x₂₁ = ODE_fase2_x₂ .- ODE_fase2_x₁

# ╔═╡ 9301f9f3-78cc-40b4-9507-cbbcba57d4f4
md"###### Velocidad relativa ``\text{v₂}\,(t) - \text{v₁}\,(t)``"

# ╔═╡ 5162a911-7352-4172-8b7e-7b515535bc79
ODE_fase2_v₂₁ = ODE_fase2_v₂ .- ODE_fase2_v₁

# ╔═╡ 00915783-c7ea-4652-b9bc-1d4518cd4d9d
md"###### Aceleración ``\text{a₁}\,(t)``"

# ╔═╡ c0692868-4b3f-4365-85d3-e6e5dbafdc45
ODE_fase2_a₁ = (k/m₁) .* (ODE_fase2_x₂₁ .- L₀)

# ╔═╡ a95f90c2-a578-4a97-9409-381cb354c605
md"###### Aceleración ``\text{a₂}\,(t)``"

# ╔═╡ c09409f3-37a1-4741-ab0f-e55cc34d0296
ODE_fase2_a₂ = -(k/m₂) .* (ODE_fase2_x₂₁ .- L₀)

# ╔═╡ 0015876b-e849-473f-a8cb-43de7304c9c4
md"###### Aceleración relativa ``\text{x₂}\,(t) - \text{x₁}\,(t)``"

# ╔═╡ 35b44d22-9f18-4f0e-aa1e-cf4b20709a01
ODE_fase2_a₂₁ = ODE_fase2_a₂ .- ODE_fase2_a₁

# ╔═╡ 5be82952-ac40-47d3-acea-54f3a3863a6e
md"###### Posición del centro de masas ``\text{Xcm}(t)``"

# ╔═╡ cc75893e-cb57-43b8-840b-f1abccab6806
ODE_fase2_Xcm = (m₁ .* ODE_fase2_x₁ .+ m₂ .* ODE_fase2_x₂)./(m₁+m₂)

# ╔═╡ 39e7f82f-1328-4f2e-8cc9-0a1214288235
md"###### Velocidad del centro de masas ``\text{Vcm}(t)``"

# ╔═╡ a14e383f-2ecf-49de-a01f-2e284d90dbc3
ODE_fase2_Vcm = (m₁ .* ODE_fase2_v₁ .+ m₂ .* ODE_fase2_v₂)./(m₁+m₂)

# ╔═╡ 7d1f9bf5-5798-44c7-b19d-4685a4695b1c
md"###### Aceleración del centro de masas ``\text{Acm}(t)``"

# ╔═╡ a171d3df-b781-4f92-89b6-43d5e3ceeaaa
ODE_fase2_Acm = (m₁.*ODE_fase2_a₁ .+ m₂.*ODE_fase2_a₂)/(m₁+m₂)

# ╔═╡ 7a46b541-1ef4-4247-aff3-707fed9032b8
md"###### Elongación del muelle ``x(t)``"

# ╔═╡ 486191ba-2176-4711-a2ac-b93203be369d
ODE_fase2_elong = L₀ .- (ODE_fase2_x₂ .- ODE_fase2_x₁)

# ╔═╡ 1ca0c73a-e697-487e-961d-f1dc5b61527e
md"###### Energía cinética de las masas ``Ec(t)``"

# ╔═╡ bf7310de-4f5e-4021-889c-a28c6249b3ac
ODE_fase2_Ec = 0.5 .* (m₁ .* ODE_fase2_v₁.^2 .+ m₂ .* ODE_fase2_v₂.^2)

# ╔═╡ 26f99d09-6000-4e6f-a3ed-5ab88b88ecd1
md"###### Energía potencial del muelle ``Ec(t)``"

# ╔═╡ 3acfcf24-2a22-466f-a9c7-accd3b68218e
ODE_fase2_Ep = 0.5 .* k .* ODE_fase2_elong.^2

# ╔═╡ e4d73407-7d74-4f1a-9f9d-4b60fd3d1028
md"###### Energía mecánica del sistema ``E(t)``"

# ╔═╡ 040a060e-aa2d-4144-aa5e-c9ee908f2e6c
ODE_fase2_E = ODE_fase2_Ec .+ ODE_fase2_Ep

# ╔═╡ e607f59c-ec4b-4d7e-bf78-8cba8ce1df50
md"#### Gráficas (fase 2)"

# ╔═╡ f140b166-65c9-444f-a8aa-d4fb6c6ea1e0
begin
	plot_fase2_x = plot(title = "Posiciones x₁, x₂, x₂-x₁, Xcm (fase 2)", 					            xaxis = "tiempo (s)", yaxis = "posición (m)",
						legend_columns = 2, legend = :outerbottom
					   );
	plot!(plot_fase2_x, ODE_fase2_t, ODE_fase2_x₁,
		  linecolor = :blue, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "x₁(t) ODE"
		 );
	plot!(plot_fase2_x, ODE_fase2_t, ODE_fase2_x₂,
		  linecolor = :red, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "x₂(t) ODE"
		 );
	plot!(plot_fase2_x, ODE_fase2_t, ODE_fase2_x₂₁,
		  linecolor = :purple, alpha=0.5, yguidefontcolor = :black,
		  linewidth = 3.0,
		  label = "x₂(t) - x₁(t) ODE"
		 );
		plot!(plot_fase2_x, ODE_fase2_t, ODE_fase2_Xcm,
		  linecolor = :grey, yguidefontcolor = :black,
		  linewidth = 2.0,
		  label = "Xcm(t) ODE"
		 )
end

# ╔═╡ dabc98e3-d5f1-47fb-85c8-d2d847e3cef3
begin
	plot_fase2_v = plot(title = "Velocidades (fase 2)",
						xaxis = "tiempo (s)", yaxis = "velocidad (m/s)",
					    legend_columns = 2, legend = :outerbottom);
	plot!(plot_fase2_v, ODE_fase2_t, ODE_fase2_v₁,
		  linecolor = :blue, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "v₁(t) ODE"
		 );
	plot!(plot_fase2_v, ODE_fase2_t, ODE_fase2_v₂,
		  linecolor = :red, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "v₂(t) ODE"
		 );
	plot!(plot_fase2_v, ODE_fase2_t, ODE_fase2_v₂₁,
		  linecolor = :purple, alpha=0.5, yguidefontcolor = :black,
		  linewidth = 3.0,
		  label = "v₂(t) - v₁(t) ODE"
		 );
	plot!(plot_fase2_v, ODE_fase2_t, ODE_fase2_Vcm,
		  linecolor = :grey, yguidefontcolor = :black,
		  linewidth = 2.0,
		  label = "Vcm(t) ODE"
		 )	
end

# ╔═╡ e68e77dc-da9a-4bee-8968-e554477bf92c
begin
	plot_fase2_a = plot(title = "Aceleraciones (fase 2)",
						xaxis = "tiempo (s)", yaxis = "aceleración (m/s2)",
					    legend_columns = 2, legend = :outerbottom);
	plot!(plot_fase2_a, ODE_fase2_t, ODE_fase2_a₁,
		  linecolor = :blue, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "a₁(t) ODE"
		 );
	plot!(plot_fase2_a, ODE_fase2_t, ODE_fase2_a₂,
		  linecolor = :red, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "a₂(t) ODE"
		 );
	plot!(plot_fase2_a, ODE_fase2_t, ODE_fase2_a₂₁,
		  linecolor = :purple, alpha=0.5, yguidefontcolor = :black,
		  linewidth = 3.0,
		  label = "a₂(t) - a₁(t) ODE"
		 );
	plot!(plot_fase2_a, ODE_fase2_t, ODE_fase2_Acm,
		  linecolor = :grey, yguidefontcolor = :black,
		  linewidth = 2.0,
		  label = "Acm(t) ODE"
		 )
end

# ╔═╡ 6cdac429-dc25-4536-8508-391c0e9b6668
md"#### Gráficas (Total)"

# ╔═╡ c000c749-125f-4936-877f-869cafad9ed5
begin
	plot_x = plot(title = "Posiciones (total)",
				  xaxis = "tiempo (s)", yaxis = "posición (m)",
				  legend_columns = 2, legend = :outerbottom);
	plot!(plot_x, ODE_fase1_t, ODE_fase1_x₁,
		  linecolor = :red, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "x₁(t) ODE fase 1"
		 );
	plot!(plot_x, ODE_fase2_t, ODE_fase2_x₁,
		  linecolor = :red, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "x₁(t) ODE fase 2"
		 );
	plot!(plot_x, ODE_fase1_t, ODE_fase1_x₂,
		  linecolor = :blue, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "x₂(t) ODE fase 1"
		 );
	plot!(plot_x, ODE_fase2_t, ODE_fase2_x₂,
		  linecolor = :blue, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "x₂(t) ODE fase 2"
		 )
end

# ╔═╡ 92d93ec3-3054-41cb-8569-857618d8f742
begin
	plot_v = plot(title = "Velocidades (total)",
				  xaxis = "tiempo (s)", yaxis = "velocidad (m/s)",
				  legend_columns = 2, legend = :outerbottom
 			);
	plot!(plot_v, ODE_fase1_t, ODE_fase1_v₁,
		  linecolor = :red, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "v₁(t) ODE fase 1"
		 );
	plot!(plot_v, ODE_fase2_t, ODE_fase2_v₁,
		  linecolor = :red, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "v₁(t) ODE fase 2"
		 );
	plot!(plot_v, ODE_fase1_t, ODE_fase1_v₂,
		  linecolor = :blue, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "v₂(t) ODE fase 2"
		 );
	plot!(plot_v, ODE_fase2_t, ODE_fase2_v₂,
		  linecolor = :blue, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "v₂(t) ODE fase 1"
		 );
end

# ╔═╡ 0ad597cf-cefc-4009-9133-46780a7258f0


# ╔═╡ 7db4a7af-b79b-4ff5-9481-a9a5fa2e71b1
begin
	plot_a = plot(title = "Aceleraciones (total)",
				  xaxis = "tiempo (s)", yaxis = "aceleración (m/s²)",
				  legend_columns = 2, legend = :outerbottom
 			);
	plot!(plot_a, ODE_fase1_t, ODE_fase1_a₁,
		  linecolor = :red, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "a₁(t) ODE fase 1"
		 );
	plot!(plot_a, ODE_fase2_t, ODE_fase2_a₁,
		  linecolor = :red, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "a₁(t) ODE fase 2"
		 );
	plot!(plot_a, ODE_fase1_t, ODE_fase1_a₂,
		  linecolor = :blue, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "v₂(t) ODE fase 1"
		 );
	plot!(plot_a, ODE_fase2_t, ODE_fase2_a₂,
		  linecolor = :blue, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "a₂(t) ODE fase 2"
		 )
end

# ╔═╡ fce3d85c-b507-43a5-a33c-457948d00f86
begin
	plot_energ = plot(title = "Energía cinética, potencial y mecánica (total)",
				  xaxis = "tiempo (s)", yaxis = "energía (J)",
				  legend_columns = 2, legend = :outerbottom
 			);
	plot!(plot_energ, ODE_fase1_t, ODE_fase1_Ec,
		  linecolor = :red, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "Ec(t) ODE fase 1"
		 );
	plot!(plot_energ, ODE_fase2_t, ODE_fase2_Ec,
		  linecolor = :red, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "Ec(t) ODE fase 2"
		 );
	plot!(plot_energ, ODE_fase1_t, ODE_fase1_Ep,
		  linecolor = :blue, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "Ep(t) ODE fase 1"
		 );
	plot!(plot_energ, ODE_fase2_t, ODE_fase2_Ep,
		  linecolor = :blue, yguidefontcolor = :black,
		  linewidth = 5.0,
		  label = "Ep(t) ODE fase 2"
		 );
	plot!(plot_energ, ODE_fase1_t, ODE_fase1_E,
		  linecolor = :black, alpha = 0.5, yguidefontcolor = :black,
		  linewidth = 2.0,
		  label = "E(t) ODE fase 1"
		 );
	plot!(plot_energ, ODE_fase2_t, ODE_fase2_E,
		  linecolor = :black, yguidefontcolor = :black,
		  linewidth = 2.0,
		  label = "E(t) ODE fase 2"
		 )
end

# ╔═╡ 7cd9047e-2b93-4ef3-9d7e-e11ef4b2196f
md"#### Vectores totales"

# ╔═╡ 9f7bdcef-6373-4338-912a-04585df5fcfd
# Vector tiempo
ODE_t = vcat(ODE_fase1_t, ODE_fase2_t);

# ╔═╡ 4ffa4437-0a12-4b8e-a970-3cb6ef136d75
begin
	# Vectores posiciones
	ODE_x₁ = vcat(ODE_fase1_x₁, ODE_fase2_x₁);
	ODE_x₂ = vcat(ODE_fase1_x₂, ODE_fase2_x₂);
end

# ╔═╡ 357df361-3f01-4fbe-bfab-897a5b8ee563
begin
	# Vectores velocidades
	ODE_v₁ = vcat(ODE_fase1_v₁, ODE_fase2_v₁);
	ODE_v₂ = vcat(ODE_fase1_v₂, ODE_fase2_v₂);
end

# ╔═╡ bdb70ea9-9824-4d32-b37e-84ad0a0a73b3
begin
	# Vectores aceleraciones
	ODE_a₁ = vcat(ODE_fase1_a₁, ODE_fase2_a₁);
	ODE_a₂ = vcat(ODE_fase1_a₂, ODE_fase2_a₂);
end

# ╔═╡ 4c5c02d3-d9fd-4f12-b1dc-2f013fe6db4e
xmin = min(minimum(ODE_x₁), minimum(ODE_x₂))

# ╔═╡ 5c5806a0-338b-4370-a0e0-e8e60fb600a0
xmax = max(maximum(ODE_x₁), maximum(ODE_x₂))

# ╔═╡ 5428d08b-4550-41df-92c1-0f1580285a76
p₁ = 100 .* (ODE_x₁ .- xmin)./(xmax - xmin)

# ╔═╡ dd68bd9c-579d-4ea6-97b2-660859fd63df
p₂ = 100 .* (ODE_x₂ .- xmin)./(xmax - xmin)

# ╔═╡ 87bd5b0f-ba77-4364-af2a-26476df836bb
# ╠═╡ disabled = true
#=╠═╡
simulacion = true
  ╠═╡ =#

# ╔═╡ 18ca6a2c-eea2-422b-963a-760ae0780c3a
# ╠═╡ disabled = true
#=╠═╡
begin
	simulacion
	plot_name = "simul"
	len = length(p₁)
	
	plot_name = plot(legend = :none,
					 #axis = ([], false),
					 size = (1024, 150),
					 xlims = (-2.0, 102.0),
					 ylims = (0.0, 100.0)
	);

#	scatter!(plot_name,
#			 [(x₁[1], 20.0), (x₂[1], 20.0)],
#			 markersize = 5.0,
#			 markercolor = :blue
#	);

	for i in 1:len
		scatter!(plot_name,
				[(p₁[i], 20.0), (p₂[i], 20.0)],
				 markersize = 5.0,
				 markercolor = :blue
		)
	end

	#simula_muelle(plot_name, p₁, p₂)
end
  ╠═╡ =#

# ╔═╡ 55baab75-4c3c-45e1-a7ee-0d7d56cd4052
# ╠═╡ disabled = true
#=╠═╡
function simula_muelle(plot_name, x₁, x₂)
	len = length(x₁)
	
	plot_name = plot(legend = :none,
					 #axis = ([], false),
					 size = (1024, 150),
					 xlims = (-2.0, 102.0),
					 ylims = (0.0, 100.0)
	);

#	scatter!(plot_name,
#			 [(x₁[1], 20.0), (x₂[1], 20.0)],
#			 markersize = 5.0,
#			 markercolor = :blue
#	);

	for i in 1:len
		scatter!(plot_name,
				[(x₁[i], 20.0), (x₂[i], 20.0)],
				 markersize = 5.0,
				 markercolor = :blue
		)
	end
end
  ╠═╡ =#

# ╔═╡ a80c8645-31cf-470c-9900-2a725a024cd6
# ╠═╡ disabled = true
#=╠═╡
function simula_muelle(plot_name, x₁, x₂)
	len = length(x₁)
	plot_name = plot(legend = :none,
					 axis = ([], false),
					 xlims = (-2.0, 102.0),
					 ylims = (0.0, 100.0)
	);
	ac
	scatter!(plot_name, [(x₁[1], 20.0), (x₂[1], 20.0)],
			 markersize = 5.0,
			 markercolor = :blue
	)
end
  ╠═╡ =#

# ╔═╡ 1858167c-c799-4c02-ae80-b17af06f0fda
# ╠═╡ disabled = true
#=╠═╡
#	for i in 2:len
#		scatter!(plot_name, [(x₁[i], 20.0), (x₂[i], 20.0)]
#		)

#		#scatter!(plot_name, [(x₁[i-1], 20.0), (x₂[i-1], 20.0)], alpha = 0.0)
#	end
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DifferentialEquations = "~7.16.1"
Plots = "~1.40.19"
PlutoUI = "~0.7.71"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "672213a6e48d05401b183bc4f298612970e46b27"

[[deps.ADTypes]]
git-tree-sha1 = "60665b326b75db6517939d0e1875850bc4a54368"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.17.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "3b86719127f50670efe356bc11073d84b4ed7a5d"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.42"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "f7817e2e585aa6d924fd714df1e2a84be7896c60"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.3.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.AlmostBlockDiagonals]]
deps = ["ConcreteStructs"]
git-tree-sha1 = "743abe5e5fe8cff96dad4123f263c0d8eee281c0"
uuid = "a95523ee-d6da-40b5-98cc-27bc505739d5"
version = "0.1.10"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "9606d7832795cbef89e06a550475be300364a8aa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.19.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "120e392af69350960b1d3b89d41dcc1d66543858"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.11.2"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "e35c672b239c5105f597963c33e740eeb46cf0ab"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.9.4"

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"
    CliqueTreesExt = "CliqueTrees"

    [deps.BandedMatrices.weakdeps]
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BoundaryValueDiffEq]]
deps = ["ADTypes", "BoundaryValueDiffEqAscher", "BoundaryValueDiffEqCore", "BoundaryValueDiffEqFIRK", "BoundaryValueDiffEqMIRK", "BoundaryValueDiffEqMIRKN", "BoundaryValueDiffEqShooting", "DiffEqBase", "FastClosures", "ForwardDiff", "LinearAlgebra", "Reexport", "SciMLBase"]
git-tree-sha1 = "d6ec33e4516b2e790a64128afdb54f3b536667a7"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "5.18.0"

    [deps.BoundaryValueDiffEq.extensions]
    BoundaryValueDiffEqODEInterfaceExt = "ODEInterface"

    [deps.BoundaryValueDiffEq.weakdeps]
    ODEInterface = "54ca160b-1b9f-5127-a996-1867f4bc2a2c"

[[deps.BoundaryValueDiffEqAscher]]
deps = ["ADTypes", "AlmostBlockDiagonals", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield"]
git-tree-sha1 = "47c833c459738a3f27c5b458ecf7832a4731ef4d"
uuid = "7227322d-7511-4e07-9247-ad6ff830280e"
version = "1.8.0"

[[deps.BoundaryValueDiffEqCore]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "ConcreteStructs", "DiffEqBase", "ForwardDiff", "LineSearch", "LinearAlgebra", "Logging", "NonlinearSolveFirstOrder", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings"]
git-tree-sha1 = "b7b4d8cc80f116eab2eb6124dba58ea7aef31b85"
uuid = "56b672f2-a5fe-4263-ab2d-da677488eb3a"
version = "1.11.1"

[[deps.BoundaryValueDiffEqFIRK]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "325e6981a414cfa5181218936c23f0e16dee8f08"
uuid = "85d9eb09-370e-4000-bb32-543851f73618"
version = "1.9.0"

[[deps.BoundaryValueDiffEqMIRK]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "da6ae5e564ad06ced4d7504929c58130558007dd"
uuid = "1a22d4ce-7765-49ea-b6f2-13c8438986a6"
version = "1.9.0"

[[deps.BoundaryValueDiffEqMIRKN]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "609c2d03ea024df0d475fee483b93cf0e87c29d6"
uuid = "9255f1d6-53bf-473e-b6bd-23f1ff009da4"
version = "1.8.0"

[[deps.BoundaryValueDiffEqShooting]]
deps = ["ADTypes", "ArrayInterface", "BandedMatrices", "BoundaryValueDiffEqCore", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastAlmostBandedMatrices", "FastClosures", "ForwardDiff", "LinearAlgebra", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "ba9bd1f31b58bfd5e48a56da0a426bcbd3462546"
uuid = "ed55bfe0-3725-4db6-871e-a1dc9f42a757"
version = "1.9.0"

[[deps.BracketingNonlinearSolve]]
deps = ["CommonSolve", "ConcreteStructs", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "a9014924595b7a2c1dd14aac516e38fa10ada656"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.3.0"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveChainRulesCoreExt = ["ChainRulesCore", "ForwardDiff"]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "a656525c8b46aa6a1c76891552ed5381bb32ae7b"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.30.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "0037835448781bb46feb39866934e243886d756a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "FastBroadcast", "ForwardDiff", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqRosenbrock", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SimpleUnPack", "SymbolicIndexingInterface"]
git-tree-sha1 = "d9b1e66070ce15bc2b9c3d5af6b94f693fc03ba4"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.58.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ConcreteStructs", "DataStructures", "DocStringExtensions", "EnumX", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "Parameters", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "1b1e070e57681d1176d99a5ec455717e24686612"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.183.2"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDistributionsExt = "Distributions"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseForwardDiffExt = ["ForwardDiff"]
    DiffEqBaseGTPSAExt = "GTPSA"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseMooncakeExt = "Mooncake"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["ConcreteStructs", "DataStructures", "DiffEqBase", "DifferentiationInterface", "Functors", "LinearAlgebra", "Markdown", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "397ef6fffcf418ba55264ba785b032b8a136903b"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.9.0"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "516d553f5deee7c55b2945b5edf05b6542837887"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.24.1"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "NonlinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "afdc7dfee475828b4f0286d63ffe66b97d7a3fa7"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.16.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "16946a4d305607c3a4af54ff35d56f0e9444ed0e"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.7"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "3e6d038b77f22791b8e3472b7c633acea1ecac06"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.120"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EnumX]]
git-tree-sha1 = "bddad79635af6aec424f53ed8aad5d7555dc6f00"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.5"

[[deps.EnzymeCore]]
git-tree-sha1 = "8272a687bca7b5c601c0c24fc0c71bff10aafdfd"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.12"
weakdeps = ["Adapt"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7bb1361afdb33c7f2b085aa49ea8fe1b0fb14e58"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.1+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "cae251c76f353e32d32d76fae2fea655eab652af"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.27.0"
weakdeps = ["StaticArrays"]

    [deps.ExponentialUtilities.extensions]
    ExponentialUtilitiesStaticArraysExt = "StaticArrays"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "83dc665d0312b41367b7263e8a4d172eac1897f4"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.4"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3a948313e7a41eb1db7a1e733e6335f17b4ab3c4"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "7.1.1+0"

[[deps.FastAlmostBandedMatrices]]
deps = ["ArrayInterface", "ArrayLayouts", "BandedMatrices", "ConcreteStructs", "LazyArrays", "LinearAlgebra", "MatrixFactorizations", "PrecompileTools", "Reexport"]
git-tree-sha1 = "9482a2b4face8ade73792c23a54796c79ed1bcbf"
uuid = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
version = "0.1.5"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "fd923962364b645f3719855c88f7074413a6ad92"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.0.2"

[[deps.FastPower]]
git-tree-sha1 = "5f7afd4b1a3969dc34d692da2ed856047325b06e"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.1.3"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerMooncakeExt = "Mooncake"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "31fd32af86234b6b71add76229d53129aa1b87a9"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.28.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "ce15956960057e9ff7f1f535400ffa14c92429a4"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.1.0"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "1828eb7275491981fa5f1752a5e126e8f26f8741"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.17"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "27299071cc29e409488ada41ec7643e0ab19091f"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.17+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "f88e0ba1f6b42121a7c1dfe93a9687d8e164c91b"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.5"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "35fbd0cefb04a516104b8e183ce0df11b70a3f1a"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.84.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7a98c6502f4632dbe9fb1973a4244eaa3324e84d"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed5e9c58612c4e081aecdb6e1a479e18462e041e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e95866623950267c1e4878846f848d94810de475"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.2+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "f5b57507a36f05509e72120aa84d5c3747dbd70e"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.17.0"

    [deps.JumpProcesses.extensions]
    JumpProcessesKernelAbstractionsExt = ["Adapt", "KernelAbstractions"]

    [deps.JumpProcesses.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "b94257a1a8737099ca40bc7271a8b374033473ed"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.1"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "52e1296ebbde0db845b356abbbe67fb82a0a116c"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.9"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "SparseArrays"]
git-tree-sha1 = "76627adb8c542c6b73f68d4bfd0aa71c9893a079"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "2.6.2"

    [deps.LazyArrays.extensions]
    LazyArraysBandedMatricesExt = "BandedMatrices"
    LazyArraysBlockArraysExt = "BlockArrays"
    LazyArraysBlockBandedMatricesExt = "BlockBandedMatrices"
    LazyArraysStaticArraysExt = "StaticArrays"

    [deps.LazyArrays.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "706dfd3c0dd56ca090e86884db6eda70fa7dd4af"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d3c8af829abaeba27181db4acb485b18d15d89c6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.1+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "97d502765cc5cf3a722120f50da03c2474efce04"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.4"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "4adee99b7262ad2a1a4bbbc59d993d24e55ea96f"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.4.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "LazyArrays", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "Setfield", "StaticArraysCore", "UnPack"]
git-tree-sha1 = "0f1a02cea457a2e26b67e105aa7ee549419c2550"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.37.0"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveBLISExt = ["blis_jll", "LAPACK_jll"]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveCUSOLVERRFExt = ["CUSOLVERRF", "SparseArrays"]
    LinearSolveCliqueTreesExt = ["CliqueTrees", "SparseArrays"]
    LinearSolveEnzymeExt = "EnzymeCore"
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    CUSOLVERRF = "a8cc9031-bad2-4722-94f5-40deabb4245c"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    LAPACK_jll = "51474c39-65e3-53ba-86ba-03b1b862ec14"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "16a726dba99685d9e94c8d0a8f655383121fc608"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "3.0.1"
weakdeps = ["BandedMatrices"]

    [deps.MatrixFactorizations.extensions]
    MatrixFactorizationsBandedMatricesExt = "BandedMatrices"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "25a6638571a902ecfb1ae2a18fc1575f86b1d4df"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.10.0"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DiffEqBase", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SparseArrays", "SparseMatrixColorings", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "d2ec18c1e4eccbb70b64be2435fc3b06fbcdc0a1"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.10.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = ["NLsolve", "LineSearches"]
    NonlinearSolvePETScExt = ["PETSc", "MPI"]
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSundialsExt = "Sundials"

    [deps.NonlinearSolve.weakdeps]
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Sundials = "c3572dad-4567-51f8-b174-8c6c989267f4"

[[deps.NonlinearSolveBase]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnzymeCore", "FastClosures", "LinearAlgebra", "Markdown", "MaybeInplace", "Preferences", "Printf", "RecursiveArrayTools", "SciMLBase", "SciMLJacobianOperators", "SciMLOperators", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "1d42a315ba627ca0027d49d0efb44e3d88db24aa"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "1.14.0"
weakdeps = ["BandedMatrices", "DiffEqBase", "ForwardDiff", "LineSearch", "LinearSolve", "SparseArrays", "SparseMatrixColorings"]

    [deps.NonlinearSolveBase.extensions]
    NonlinearSolveBaseBandedMatricesExt = "BandedMatrices"
    NonlinearSolveBaseDiffEqBaseExt = "DiffEqBase"
    NonlinearSolveBaseForwardDiffExt = "ForwardDiff"
    NonlinearSolveBaseLineSearchExt = "LineSearch"
    NonlinearSolveBaseLinearSolveExt = "LinearSolve"
    NonlinearSolveBaseSparseArraysExt = "SparseArrays"
    NonlinearSolveBaseSparseMatrixColoringsExt = "SparseMatrixColorings"

[[deps.NonlinearSolveFirstOrder]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "3f1198ae5cbf21e84b8251a9e62fa1f888f3e4cb"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.7.0"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "DiffEqBase", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "40dfaf1bf74f1f700f81d0002d4dd90999598eb2"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.8.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "DiffEqBase", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "84de5a469e119eb2c22ae07c543dc4e7f7001ee7"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.3.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveSpectralMethods.extensions]
    NonlinearSolveSpectralMethodsForwardDiffExt = "ForwardDiff"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "f1a7e086c677df53e064e0fdd2c9d0b0833e3f6e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.5.0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2ae7d4ddec2e13ad3bddf5c0796f7547cf682391"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.2+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "EnumX", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "61942645c38dd2b5b78e2082c9b51ab315315d10"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.13.2"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c392fc5dd032381919e3b22dd32d6443760ce7ea"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.5.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "55c21fdb4626037cdbcb04fec3afa192345a24de"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.101.0"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "09aae1486c767caa6bce9de892455cbdf5a6fbc8"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.5.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "ce8db53fd1e4e41c020fd53961e7314f75e4c21c"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.10.1"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "e579c9a4f9102e82da3d97c349a74d6bc11cf8dc"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.30.0"

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"
    OrdinaryDiffEqCoreMooncakeExt = "Mooncake"

    [deps.OrdinaryDiffEqCore.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.OrdinaryDiffEqDefault]]
deps = ["ADTypes", "DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport", "SciMLBase"]
git-tree-sha1 = "7d5ddeee97e1bdcc848f1397cbc3d03bd57f33e7"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.8.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SciMLOperators", "SparseArrays", "SparseMatrixColorings", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "4c270747152db513dd80bfd5f2f9df48befff28a"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.14.0"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "4c0633f587395d7aaec0679dc649eb03fcc74e73"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.4.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["ADTypes", "DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "3b81416ff11e55ea0ae7b449efc818256d9d450b"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.8.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "ee2cba2533e9faf71b09a319a910d4886931e7a6"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.8.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastGaussQuadrature", "FastPower", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "b968d66de3de5ffcf18544bc202ca792bad20710"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.16.0"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "815b54211201ec42b8829e0275ab3c9632d16cbe"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.4.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "fe750e4b8c1b1b9e1c1319ff2e052e83ad57b3ac"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.5.0"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "42096f72136078fa02804515f1748ddeb1f0d47d"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.5.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "a5dcd75959dada0005b1707a5ca9359faa1734ba"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.7.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "925fc0136e8128fd19abf126e9358ec1f997390f"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.6.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "3cc4987c8e4725276b55a52e08b56ded4862917e"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.6.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "9291cdfd2e8c91e900c48d71d76618de47daeede"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.6.0"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "b05226afc8fa6b8fc6f2258a89987b4f5bd0db4e"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.14.1"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "c90aa7fa0d725472c4098096adf6a08266c2f682"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.4.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "9d599d2eafdf74ab26ea6bf3feb28183a2ade143"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.6.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport", "SciMLBase"]
git-tree-sha1 = "8e35132689133255be6d63df4190b5fc97b6cf2b"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.4.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "63fb643a956b27cd0e33a3c6d910c3c118082e0f"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.4.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "a31c41f9dbea7c7179c6e544c25c7e144d63868c"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.5.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "a06c1263d71ea42a1881b4d49c8a087035d4a3ff"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.16.1"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "20caa72c004414435fb5769fadb711e96ed5bcd4"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.7.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "af955c61407631d281dd4c2e8331cdfea1af49be"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.6.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqStabilizedRK", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "75abe7462f4b0b2a2463bb512c8a5458bbd39185"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.6.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "7e94d3d1b3528b4bcf9e0248198ee0a2fd65a697"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.4.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "e8dd5ab225287947016dc144a5ded1fb83885638"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.7.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "778c7d379265f17f40dbe9aaa6f6a2a08bc7fa3e"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.5.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "185578fa7c38119d4318326f9375f1cba0f0ce53"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "f07c06228a1c670ae4c87d1276b92c7c597fdda0"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.35"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "275a9a6d85dc86c24d03d1837a0010226a96f540"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.3+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "0c5a5b7e440c008fe31416a3ac9e0d2057c81106"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.19"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

[[deps.PoissonRandom]]
deps = ["LogExpFunctions", "Random"]
git-tree-sha1 = "c1ea45aa9f209fe97192afa233907bc4e551c8aa"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.6"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6f7cd22a802094d239824c57d94c8e2d0f7cfc7d"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.18"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "9b4ee15d1fc68654031964a7af0c914c898e35a7"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.33"

    [deps.PreallocationTools.extensions]
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "eb38d376097f47316fe089fc62cb7c6d85383a52"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "e1d5e16d0f65762396f9ca4644a5f4ddab8d452b"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "dbe5fd0b334694e905cb9fda73cd8554333c46e2"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.1"

[[deps.RandomNumbers]]
deps = ["Random"]
git-tree-sha1 = "c6ec94d2aaba1ab2ff983052cf6a606ca5985902"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.6.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "96bef5b9ac123fff1b379acf0303cf914aaabdfd"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.37.1"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "86a8a8b783481e1ea6b9c91dd949cb32191f8ab4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.15"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "6c22e097cdc0ea9f64fa1fd9b9ba90bc208264ee"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.113.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReverseDiffExt = "ReverseDiff"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "3414071e3458f3065de7fa5aed55283b236b4907"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.8"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "78ac1b947205b07973321f67f17df8fbe6154ac9"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.6.0"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLStructures]]
deps = ["ArrayInterface"]
git-tree-sha1 = "566c4ed301ccb2a44cbd5a27da5f885e0ed1d5df"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "09d986e27a606f172c5b6cffbd8b8b2f10bf1c75"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.7.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveDiffEqBaseExt = "DiffEqBase"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffEqBase = "2b5f629d-d688-5b77-993f-72d75c75574e"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "339efef69fda0cccf14c06a483561527e9169b8f"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "1.0.1"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "9de43e0b9b976f1019bf7a879a686c4514520078"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.21"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = "CUDA"
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "41852b8679f78c8d8961eeadc8f62cef861a52e3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "95af145932c2ed859b63329952ce8d633719f091"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.3"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools"]
git-tree-sha1 = "f737d444cb0ad07e61b3c1bef8eb91203c321eff"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.2.0"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

    [deps.StaticArrayInterface.weakdeps]
    OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "cbea8a6bd7bed51b1619658dec70035e07b8502f"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.14"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9d72a13a3f4dd3795a195ac5a44d7d6ff5f552ff"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.1"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2c962245732371acd51700dbb268af311bddd719"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.6"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "8e45cecc66f3b42633b8ce14d431e8e57a3e242e"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.SteadyStateDiffEq]]
deps = ["ConcreteStructs", "DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NonlinearSolveBase", "Reexport", "SciMLBase"]
git-tree-sha1 = "66a028f9a2bb44d0f6de0814a2b9840af548143a"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "2.5.0"

[[deps.StochasticDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FastPower", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SparseArrays", "StaticArrays", "UnPack"]
git-tree-sha1 = "c3a55a2a1e180e249a0550d30a58c700487aa7ef"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.81.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "83151ba8065a73f53ca2ae98bc7274d817aa30f2"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.8"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "PrecompileTools", "Reexport", "SciMLBase", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "7c7a7ee705724b3c80d5451ac49779db36c6f758"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.28.0"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "SuiteSparse_jll", "libblastrampoline_jll"]
git-tree-sha1 = "91db7ed92c66f81435fe880947171f1212936b14"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.3+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "93104ca226670c0cb92ba8bc6998852ad55a2d4c"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.43"

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "6258d453843c466d84c17a58732dda5deeb8d3af"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.24.0"
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "af305cc62419f9bd61b6644d19170a4d258c7967"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.7.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fee71455b0aaa3440dfdd54a9a36ccef829be7d4"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.1+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "9caba99d38404b285db8801d5c45ef4f4f425a6d"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.1+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "c5bf2dad6a03dfef57ea0a170a1fe493601603f2"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.5+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4bba74fa59ab0755167ad24f98800fe5d727175b"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.12.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "07b6a107d926093898e82b3b1db657ebe33134ec"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.50+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "fbf139bce07a534df0e699dbb5f5cc9346f95cc1"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.9.2+0"
"""

# ╔═╡ Cell order:
# ╟─e5209b22-863c-11f0-2e47-e157a5f44a52
# ╟─29d4fb45-0f1b-424f-8de6-ed277df6931f
# ╟─ad8ad00a-0f02-40bd-9329-e0f76c92bb92
# ╟─88a02306-6b34-48df-ad4f-0f232b5ec826
# ╟─29fa4220-9626-4a87-bd99-49ac4ed74df4
# ╟─8df61bb0-70ee-4744-b0b9-8a71ac519d10
# ╟─bacf9614-5d5c-4d6b-a2b3-3bb4c7dd8e9f
# ╟─cc7fe116-32e5-45da-b0c7-0fff78039593
# ╟─96c1c896-4636-4a37-a6b6-40c4c1c6e729
# ╟─b8943312-9df6-4fdd-a94e-31ccff232b63
# ╟─aaaf469c-42e3-4d23-b57b-8793db22dbb0
# ╟─7e519bb7-dc56-45df-9f44-e8d530b9a049
# ╟─21606d10-ebe1-419d-b463-7dd2cf20fe71
# ╟─930be4ad-2eee-47c7-90d5-e5ad2ccb4608
# ╟─e482e2e5-7e90-4648-9d23-3c2b97ee8adb
# ╟─83d12ae9-3c71-4210-80f9-39e592d7dfe5
# ╟─83aa9ded-a4c4-47da-99db-053c67ee587b
# ╟─f02e84e1-7ddb-4b46-b7ac-185ade9d9d83
# ╟─6756649b-92e3-44e9-ab63-05db279b45df
# ╟─87744172-a35c-4360-a016-8f55ea9a67e3
# ╟─50f66526-988c-41e1-93c4-b8899f90f3e0
# ╟─07db6630-e83e-4f14-8fce-65fdb00c2153
# ╟─bb104dc8-fd20-45fe-b0b3-044038da1cd3
# ╟─f446b922-4640-44f8-acba-388089d7331a
# ╟─18d770e4-0aca-4370-aceb-2f995aa64bc9
# ╟─0691719a-f48e-4c62-b1ac-47d906bdb5d1
# ╟─ce31b020-9e41-4885-8085-48c479556735
# ╟─f7d6c8f3-aa8e-4637-a2ee-7b0512465119
# ╟─4581a34a-225a-4960-bf84-efbf63bf0039
# ╟─fd0bada9-6c77-4725-9835-99d64db36302
# ╟─54c7087b-e18b-41ff-9ed1-cfb231bf8dbc
# ╟─14969615-02d6-48be-9d42-30ba0e2175cb
# ╟─76b1fe59-534d-4d45-86b0-f758f2a34ed9
# ╟─e937e066-2b5c-4585-99cc-987e72fd5f8f
# ╟─f21e13e0-6a37-4c69-8cf3-f3bdf6da0792
# ╟─56710ca1-077e-449a-b086-dc1e69d7f558
# ╟─d7725e84-5f41-4a15-8d0c-a79208d81023
# ╟─05a9576b-dc56-4621-a00a-9ab1033b8344
# ╟─2d7d3def-2cf1-4350-aeea-2c2b1db507cc
# ╟─7a884cb3-5f42-4c34-a5b8-81ff7553e495
# ╟─c1797c7a-1b0a-436f-b7fd-841d9b5e3bc5
# ╟─7e6e1f69-fff3-495d-b73e-43e551d9b9a6
# ╟─4fea284e-df13-4110-b587-25c2ab7c413e
# ╟─601c518b-b860-438a-8777-f370b1057488
# ╟─af858578-9c8e-4b0f-ba80-6ddcc030e840
# ╟─34b7d3fc-5467-4ec4-bc18-33f69f8779b0
# ╟─56311590-1ce5-4316-92f4-6872cbcab85b
# ╟─7da08f46-7865-47b1-86d0-07e74e32ed6f
# ╟─554258ce-0ccb-47bc-a629-9a2811cfbf6f
# ╟─29f24328-0805-4360-9222-d859987dcfb8
# ╟─891335e5-a97d-42a2-a9ee-cde6bd183f9c
# ╟─606f1855-6040-4f5d-b5b1-3db3bed2c658
# ╟─f74a28f8-cb6a-4fd9-8d6c-0cce6eb72611
# ╟─b93586da-4555-4b04-8685-bfdbeebb6da2
# ╟─bfee0f02-0c8e-463b-b5d7-f025a224acff
# ╟─2ae54de3-9506-49ac-b238-5d944089b7b7
# ╟─5f523482-86a2-4f36-9f8e-ef3a4d72e53a
# ╟─ebe29b79-59c2-4cfa-b964-9d978345bd91
# ╟─ea584508-2b28-4942-80a4-1fccd9f4e954
# ╟─a22f3e68-eed3-4182-b169-046464a87172
# ╟─7874884f-2cb1-4c71-a757-4ae5cff6b5c0
# ╟─e912655a-f2bd-4d7f-aaea-2612c25a5b3e
# ╟─b6319630-97cf-43d2-bb73-aedeeaa903a1
# ╟─dc48a744-10d9-4885-982f-ebca319f4aee
# ╟─3668c039-e55e-4bf3-959c-1613f9b8b481
# ╟─afe05bcd-357d-4869-83fd-89a398504539
# ╟─2a30b827-a4ac-4f15-914f-b2941c6ac10a
# ╟─ec6511a6-baa1-4959-8f45-4f1deb9d49fd
# ╟─958c56e3-26ef-4cf2-b3e8-ce7ff5b9d7e6
# ╟─ac175b42-b5a9-4462-805b-1de215710578
# ╟─6ec06427-b723-4fc0-8bec-8a92b74417c0
# ╟─0a3d31d0-5c25-4d18-b176-22f1eee61cd5
# ╟─95b48e57-9b22-442d-9f79-24e3ab148bf3
# ╟─f8e14404-bdc8-4f48-b278-fbcf6a4566dd
# ╟─5bcb6eec-59b9-4096-a95b-c136281baad1
# ╟─3a788f13-59b8-48a3-ab9f-48b96d2a3486
# ╟─be85ce71-110a-4dd5-bbfd-52f19784d2e0
# ╟─381e74e7-4754-4ad7-8625-520c3674c666
# ╟─1ba5ad0c-be87-4629-9c19-f81d54cad92b
# ╟─064cd216-775e-4995-b13c-ed00c8e56b51
# ╟─c296ceca-501c-4cdb-bc76-6a9f796d38cb
# ╟─e7424c03-0db4-447e-8f17-687a6f01347c
# ╟─c159a02c-2384-4c11-a1e5-03438de7cda6
# ╟─75494f1e-a3f7-4790-b364-f58fa8910e42
# ╟─299f987d-fa02-4fe7-901d-d2c40d8e785a
# ╟─88630861-b840-44cf-8230-39436b2888fc
# ╟─2736abc6-c8fa-42f9-a7f1-8c4797dea825
# ╟─2dc09545-4c46-4833-8d76-45e6383fad7d
# ╟─173325de-9ac5-48a1-912c-3806212a1796
# ╟─052bf6ad-9fcc-4037-9115-dcbb2f19d99b
# ╟─84b4fc0b-ee8d-43d0-acb9-e1b3d8146224
# ╟─17f19d51-a5a3-497f-ba26-68cdc7ff9834
# ╟─8fbdd5bb-ba8c-4ba9-a90e-ae804c50d61b
# ╟─fb43062f-fed6-4120-bc4b-2af832c16387
# ╟─40b922a7-5155-4cf6-a50c-4597a32c2bf4
# ╟─26dffbe9-d8c4-402b-9913-7172c0458b8e
# ╟─0c4207fd-3b47-4042-b821-ca24917145ed
# ╟─9c9f9325-e5e1-4ebb-9402-fdbf49de8e32
# ╟─7f4dee6a-7e72-4531-9c42-d02fd8d3bf12
# ╟─924f3ca7-a169-4297-846c-e88607d6fd54
# ╠═e89a98ac-917b-4094-871d-6f75a344567b
# ╟─28c31f91-e55f-4d3b-b27d-a0ea722a7447
# ╠═28a797d4-7e36-4283-9daa-aefae713c4a4
# ╟─5a5bcb8d-f0fb-406f-acd4-680edc09c1df
# ╟─3dccada0-19c5-4dad-9616-f73a9d00d86d
# ╟─e04a6c73-e3cf-48ab-908f-fd16798062a0
# ╟─5486689c-ab65-439f-98ce-e358e6ba24f2
# ╟─0c8af1d1-381d-4166-925b-246ba9e72d08
# ╟─f12c8235-b240-494d-91ee-0b264ba9a4c0
# ╟─7807081f-c329-4e26-b4d2-0d2848d2bc0e
# ╟─7d3ba095-977d-4108-857f-039dc50b022c
# ╟─b0534f3d-6481-4e3d-aa48-0888eeceab23
# ╟─b0142e00-5f55-4e5f-bd49-962b862f6312
# ╟─236e9c09-571e-47a9-92c1-fa9ac4e6965e
# ╟─21eaf4af-2f31-4f2a-8cbf-c3719c1db356
# ╟─d952b13c-ef40-4f46-bfcd-883f951faf74
# ╠═c3ab23be-07f4-4429-834a-a612ae221a4a
# ╟─c677f634-3133-44ce-bc5e-04c2d516f742
# ╟─5e3d0230-26be-42bb-a057-4784be1263f5
# ╟─065d6b1c-67de-4070-9c1e-774553b49168
# ╠═e0e48f7a-dcd7-4fb6-ae6a-1f3635618cfc
# ╟─e660648a-8a93-4efc-9e58-5207307263be
# ╟─18fc5377-34ae-4782-bf7e-62b741286044
# ╠═57432286-84f8-4c66-84ba-2c03c5769d5e
# ╟─83deae0a-a7ce-400c-8dbe-07bce3d8e8f8
# ╠═ebb2888e-4561-44e3-998d-60241bca7c2d
# ╟─8e45f531-a6ed-48c7-a21c-c218dae6a9f3
# ╟─433dd48b-cd87-48ef-86d6-8e32959c2fbf
# ╟─a31afff7-a77f-4b36-b576-7ef3f9f33380
# ╟─b97fb2e2-9f61-4555-b50f-94c5c69c5b71
# ╠═c2f0923e-58b9-4b2f-8c8d-f6174e6e4523
# ╟─0fb3aebe-affb-4849-a987-5fe11c900fd9
# ╠═d6eacbd0-7e55-4c8d-829f-54a4f47051b9
# ╟─f475a8f0-d55b-48aa-b5e9-fcfe358add4e
# ╠═8d131541-edc6-4cd4-b35a-4bb4981dbee6
# ╟─03b0c959-9d0d-46e4-96c2-972a0771922f
# ╠═6a887098-abd4-41ba-9e50-de4b7eae2fc9
# ╟─369d58f4-6c8d-48af-b172-751e89f4ceea
# ╠═d7ab580c-061f-4369-babf-1c7608d24073
# ╟─f76e194d-65da-48bc-ac2c-c213473f79e6
# ╠═0b53878d-9a5d-4b96-a2eb-5a2bdc1e4027
# ╟─1cfc9bbd-05cd-4d31-90e9-f8ba2858fb23
# ╠═073cd0dc-4932-455b-9f88-a350e2e9695a
# ╟─8b369410-4589-4878-996d-9cf034e1281f
# ╠═1552fbf4-a1c3-4daa-9ead-387c11dfa23a
# ╟─372f8d9a-7cb1-4ff6-8cbe-b81c57872feb
# ╠═7b531ed8-44b4-4a60-b182-24b1a0f5dc01
# ╟─80b47e79-e777-4d6c-9ef7-ffdea8fcf1ee
# ╠═0d733fd6-27cd-4a54-bced-9f579a7c3ef3
# ╟─9175746e-286f-4856-ad16-f64d479f7842
# ╠═ff64e888-2731-409c-86f9-2d10418375a7
# ╟─db21bff9-51ad-490e-ae9f-cc2deea8dfeb
# ╠═61cf579c-3807-4234-8a40-99ccd1cd54f1
# ╟─a89390ae-cbdd-4144-8f0c-ed30e8c8f8d7
# ╠═f3b6b1a2-b267-4035-9267-b8830318de17
# ╟─68cd7717-d64f-4298-9ff1-1f4ed61b5b86
# ╠═d1054a56-c567-4f82-aaac-5218e09e71cd
# ╟─5dd36860-bd12-41e3-89a2-0fec5b2ad710
# ╠═729d72ab-51f7-410f-852a-3d85cd39ff4c
# ╟─e556f528-b5b7-437b-b325-59af944e6697
# ╠═89cf5e1f-d657-45ac-82ac-edb72c1aeb3b
# ╟─b35f2726-e632-466a-8f5e-e514c4c53697
# ╠═92899b96-f19e-446f-96fa-1fc2c2d6e04d
# ╟─0cd13e27-cffa-43c0-b5de-cbc0158845e8
# ╠═244fa152-13ab-45f5-b7f4-a66067fd53d0
# ╟─ac0ca9aa-32d9-4e71-a3cf-a1f95e253a51
# ╠═13c7f234-0abf-4886-8989-48e8f2ce4b5c
# ╟─cb0ffc15-5abd-4bfa-ae0f-c05769f2f386
# ╟─2a42e1e2-5123-4888-8281-164d44229425
# ╠═b9b8faee-e7f6-4bcf-b4a7-e6cefe5e2341
# ╠═44d5951c-9562-4567-9dea-bd6425ba4676
# ╠═60db7902-2381-43e1-9a16-8d8ff33ce556
# ╠═0f89d7bf-75e9-485a-bd48-27358295c0b0
# ╟─f43f4f52-d53e-4b46-98a1-fb80cc4caa83
# ╟─2ffbe9cc-90fe-4c80-a1ca-94accea6f63a
# ╟─b622ca7d-2121-4617-bc2c-05c51b465a65
# ╠═171bf0fd-a545-4161-ae77-f3a23aa19233
# ╟─f021828f-e7d6-46a1-90cf-017345536bdd
# ╟─5b712a32-13e7-4504-bb8a-90e92d6d246a
# ╠═703907bb-a878-4823-87d1-ed00fb0ff0ed
# ╟─5c60e568-d358-42b6-aaa6-97fd1f51302a
# ╠═fbe27e04-9ff7-4b63-b91d-08d7513cce32
# ╟─84b89a14-4b20-4525-8887-c4ad35e6d74a
# ╟─feb81d34-8870-4fd0-ae46-0704d7917a9d
# ╟─2859858a-9918-450b-8851-a5995645ba84
# ╟─82f1437d-704a-4566-b70e-1061cbbf8b7a
# ╟─3ebf307c-089e-4090-b875-d599902c63f7
# ╟─2793e637-60f4-46c0-8de5-a699cf4554d0
# ╠═a73a864c-2603-4952-818c-3f48484b40cf
# ╟─0d68a866-5182-489b-bb0e-bf71e63bdb7e
# ╠═169b2496-f2f6-4966-b971-b6615ea8e241
# ╟─b6f47d8f-c1c9-4a89-a665-1157daf728f7
# ╠═6377237c-44cf-41f7-b8b8-b901ab4a1dd6
# ╟─a2448f50-744e-443c-8ee1-d2730f7b9161
# ╠═849c2bce-edb3-4796-9414-c54625aa0977
# ╟─e107ca77-2312-4fcd-a22a-58d524ac8aec
# ╠═8ca72f78-588d-4c69-a8e8-4fa16fb4b93f
# ╠═ab04a12a-b375-4038-ad6c-879033ffb3b6
# ╟─4fa8e88a-f332-4446-88f2-90abc693723f
# ╠═60f893c6-429a-4bbd-a534-b0d988afb5e4
# ╟─9301f9f3-78cc-40b4-9507-cbbcba57d4f4
# ╠═5162a911-7352-4172-8b7e-7b515535bc79
# ╟─00915783-c7ea-4652-b9bc-1d4518cd4d9d
# ╠═c0692868-4b3f-4365-85d3-e6e5dbafdc45
# ╟─a95f90c2-a578-4a97-9409-381cb354c605
# ╠═c09409f3-37a1-4741-ab0f-e55cc34d0296
# ╟─0015876b-e849-473f-a8cb-43de7304c9c4
# ╠═35b44d22-9f18-4f0e-aa1e-cf4b20709a01
# ╟─5be82952-ac40-47d3-acea-54f3a3863a6e
# ╠═cc75893e-cb57-43b8-840b-f1abccab6806
# ╟─39e7f82f-1328-4f2e-8cc9-0a1214288235
# ╠═a14e383f-2ecf-49de-a01f-2e284d90dbc3
# ╟─7d1f9bf5-5798-44c7-b19d-4685a4695b1c
# ╠═a171d3df-b781-4f92-89b6-43d5e3ceeaaa
# ╟─7a46b541-1ef4-4247-aff3-707fed9032b8
# ╠═486191ba-2176-4711-a2ac-b93203be369d
# ╠═1ca0c73a-e697-487e-961d-f1dc5b61527e
# ╠═bf7310de-4f5e-4021-889c-a28c6249b3ac
# ╟─26f99d09-6000-4e6f-a3ed-5ab88b88ecd1
# ╠═3acfcf24-2a22-466f-a9c7-accd3b68218e
# ╟─e4d73407-7d74-4f1a-9f9d-4b60fd3d1028
# ╠═040a060e-aa2d-4144-aa5e-c9ee908f2e6c
# ╟─e607f59c-ec4b-4d7e-bf78-8cba8ce1df50
# ╠═f140b166-65c9-444f-a8aa-d4fb6c6ea1e0
# ╠═dabc98e3-d5f1-47fb-85c8-d2d847e3cef3
# ╠═e68e77dc-da9a-4bee-8968-e554477bf92c
# ╟─6cdac429-dc25-4536-8508-391c0e9b6668
# ╠═c000c749-125f-4936-877f-869cafad9ed5
# ╠═92d93ec3-3054-41cb-8569-857618d8f742
# ╠═0ad597cf-cefc-4009-9133-46780a7258f0
# ╠═7db4a7af-b79b-4ff5-9481-a9a5fa2e71b1
# ╠═fce3d85c-b507-43a5-a33c-457948d00f86
# ╟─7cd9047e-2b93-4ef3-9d7e-e11ef4b2196f
# ╠═9f7bdcef-6373-4338-912a-04585df5fcfd
# ╠═4ffa4437-0a12-4b8e-a970-3cb6ef136d75
# ╠═357df361-3f01-4fbe-bfab-897a5b8ee563
# ╠═bdb70ea9-9824-4d32-b37e-84ad0a0a73b3
# ╠═4c5c02d3-d9fd-4f12-b1dc-2f013fe6db4e
# ╠═5c5806a0-338b-4370-a0e0-e8e60fb600a0
# ╠═5428d08b-4550-41df-92c1-0f1580285a76
# ╠═dd68bd9c-579d-4ea6-97b2-660859fd63df
# ╠═87bd5b0f-ba77-4364-af2a-26476df836bb
# ╠═18ca6a2c-eea2-422b-963a-760ae0780c3a
# ╠═55baab75-4c3c-45e1-a7ee-0d7d56cd4052
# ╠═a80c8645-31cf-470c-9900-2a725a024cd6
# ╠═1858167c-c799-4c02-ae80-b17af06f0fda
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
