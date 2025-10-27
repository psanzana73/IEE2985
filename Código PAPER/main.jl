include("modelo/unit_commitment.jl")
include("modelo/parametros.jl")
include("modelo/variables.jl")
include("modelo/restricciones.jl")

using .ModeloUC
modelo, par, set, var = ModeloUC.construir_modelo()
status = ModeloUC.solve_modelo(modelo)