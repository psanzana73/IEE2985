module ModeloUC

using JuMP, MathOptInterface, Gurobi
const MOI = MathOptInterface

import ..Parametros
import ..Variables
import ..Restricciones
import ..Conjuntos
import ..Objetivo

export construir_modelo, solve_modelo

function construir_modelo(data_path::AbstractString)
    
    # 1. Leer datos dinámicamente desde la ruta proporcionada por el Batch Run
    # Se capturan los 6 retornos de read_input_data
    gens, buses, branches, ibgs, scc, freq = Parametros.read_input_data(data_path)
    
    # 2. Empaquetar en la estructura 'par'
    par = (generators = gens, 
           buses      = buses, 
           impedances = branches, 
           ibgs       = ibgs, 
           scc        = scc, 
           freq       = freq)

    # 3. Inicializar el optimizador
    modelo = Model(Gurobi.Optimizer)

    # 4. Construir componentes del modelo
    set = Conjuntos.definir_conjuntos(par)
    var = Variables.definir_variables(modelo, set)
    
    Objetivo.funcion_objetivo(modelo, par, set, var)
    
    # Generar restricciones (incluyendo SCC y Frecuencia si están activas)
    Restricciones.generar_restricciones(modelo, par, set, var)

    return modelo, par, set, var
end

function solve_modelo(modelo::JuMP.Model)
    optimize!(modelo)
    
    status = MOI.get(modelo, MOI.TerminationStatus())
    
    # Opcional: Imprimir objetivo si es óptimo para feedback inmediato en consola
    if status == MOI.OPTIMAL
        println("   -> Solución Óptima encontrada. Costo = ", objective_value(modelo))
    end
    
    return status
end

end
