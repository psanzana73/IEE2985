module ModeloUC

using JuMP, MathOptInterface, Gurobi
const MOI = MathOptInterface

# Importa todos los módulos de tu modelo
import ..Parametros
import ..Conjuntos
import ..Variables
import ..Objetivo
import ..Restricciones

# Define la ruta a tus datos
const DATA_PATH = joinpath(@__DIR__, "data", "input", "datos_finales30_pen27.xlsx")

export construir_modelo, solve_modelo

# --- ELIMINADAS FUNCIONES REDUNDANTES ---
# Las funciones _bus_key y _build_Y0 se eliminaron 
# porque ya están implementadas (correctamente) 
# dentro de tu módulo Restricciones.jl

function construir_modelo(data_path::AbstractString = DATA_PATH)
    
    # CORREGIDO: Capturar los 6 parámetros, incluyendo 'freq'
    gens, buses, branches, ibgs, scc, freq = Parametros.read_input_data(data_path)
    
    # CORREGIDO: Añadir 'freq' al struct 'par'
    par = (generators = gens, 
           buses = buses, 
           impedances = branches, 
           ibgs = ibgs, 
           scc = scc, 
           freq = freq)

    modelo = Model(Gurobi.Optimizer)

    set = Conjuntos.definir_conjuntos(par)
    var = Variables.definir_variables(modelo, set)
    Objetivo.funcion_objetivo(modelo, par, set, var)
    
    # Esta función ahora añadirá las restricciones de UC y SCC
    Restricciones.generar_restricciones(modelo, par, set, var) 

    # NOTA: Tu módulo `Restricciones` actual AÚN NO CONTIENE
    # las Restricciones de Frecuencia (FreqParams). 
    # Deberás añadirlas dentro de `generar_restricciones` 
    # para que el modelo esté completo según el paper.

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
