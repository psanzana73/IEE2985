include(joinpath(@__DIR__, "parametros.jl"))
include(joinpath(@__DIR__, "variables.jl"))
include(joinpath(@__DIR__, "restricciones.jl"))

module ModeloUC

using JuMP
using MathOptInterface
using Gurobi
const MOI = MathOptInterface

import ..Parametros
import ..Conjuntos
import ..Variables
import ..Objetivo
import ..Restricciones

const DATA_PATH = joinpath(@__DIR__, "data", "input", "datos_finales30.xlsx")

export construir_modelo, solve_modelo

function construir_modelo(data_path::AbstractString = DATA_PATH)
    gens, buses, branches = Parametros.read_input_data(data_path)
    par = (generators = gens, buses = buses, impedances = branches)
    modelo = Model(Gurobi.Optimizer)

    set = Conjuntos.definir_conjuntos(par)
    var = Variables.definir_variables(modelo, set)
    Objetivo.funcion_objetivo(modelo, par, set, var)
    Restricciones.generar_restricciones(modelo, par, set, var)

    return modelo, par, set, var
end

function solve_modelo(modelo::JuMP.Model)
    optimize!(modelo)
    status = MOI.get(modelo, MOI.TerminationStatus())
    println("Estado del solver: ", status)
    if status == MOI.OPTIMAL
        println("Objetivo = ", objective_value(modelo))
    end
    return status
end

end


