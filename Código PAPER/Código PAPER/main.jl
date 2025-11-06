include("modelo/parametros.jl")
include("modelo/variables.jl")
include("modelo/restricciones.jl")

# 2. Incluir el módulo principal que los utiliza
# (Este archivo debe ir al final, ya que usa los de arriba)
include("modelo/unit_commitment.jl") 

# 3. Usar el módulo principal para construir y resolver
using .ModeloUC

println("Construyendo el modelo...")
modelo, par, set, var = ModeloUC.construir_modelo()
println("Modelo construido. Resolviendo...")

status = ModeloUC.solve_modelo(modelo)

println("Fin de la ejecución.")