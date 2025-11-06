module Conjuntos
using JuMP
export definir_conjuntos

function definir_conjuntos(par; j::Integer=1)
    G = 1:length(par.generators)
    B = 1:length(par.buses)
    C = 1:length(par.ibgs) # Conjunto de generadores IBG
    T = hasproperty(par, :T) ? par.T :
        (hasproperty(par, :time_labels) ? length(par.time_labels) : 24)

    return (GeneratorSet = G,
            BusSet       = B,
            IBGSet       = C,
            TimeSet      = 1:T,
            PiecewiseSet = 1:j,
            T            = T)
end

end

module Variables
using JuMP
export definir_variables

function definir_variables(modelo::JuMP.Model, set)
    u     = @variable(modelo, u[set.GeneratorSet,set.TimeSet], Bin) # 1 si está encendida (Commitment)
    p     = @variable(modelo, p[set.GeneratorSet,set.TimeSet] >= 0) # Potencia generada
    w     = @variable(modelo, w[set.GeneratorSet,set.TimeSet], Bin) # 1 si arranca (Startup)
    v     = @variable(modelo, v[set.GeneratorSet,set.TimeSet], Bin) # 1 si se apaga (Shutdown)
    theta = @variable(modelo, theta[set.BusSet,set.TimeSet])
    
    # variables SCC
    Z     = @variable(modelo, Z[set.BusSet, set.BusSet, set.TimeSet])  
    mu    = @variable(modelo, mu[set.BusSet, set.GeneratorSet, set.TimeSet])
    
    # Nota: 'w' es 'y' en Restricciones, 'v' es 'z'. El alias es correcto.
    return (u=u, p=p, w=w, v=v, theta=theta, Z=Z, mu=mu)
end

end

module Objetivo
using JuMP
export funcion_objetivo

function funcion_objetivo(modelo::JuMP.Model, par, set, var)
    @objective(modelo, Min,
        sum(par.generators[g].VariableCost  * var.p[g,t] +  # 1. Costo Marginal (por MW)
             par.generators[g].NoLoadCost    * var.u[g,t] +  # 2. Costo en Vacío (por hora encendida) - FALTABA
             par.generators[g].StartUpCost   * var.w[g,t] +  # 3. Costo de Arranque (por arranque) - CORREGIDO
             par.generators[g].ShutDownCost  * var.v[g,t]   # 4. Costo de Parada (por parada)
             for g in set.GeneratorSet, t in set.TimeSet)
    )
    return nothing
end

end


