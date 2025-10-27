module Conjuntos
using JuMP
export definir_conjuntos

function definir_conjuntos(par; j::Integer=1)
    G = 1:length(par.generators)
    B = 1:length(par.buses)
    T = hasproperty(par, :T) ? par.T :
        (hasproperty(par, :time_labels) ? length(par.time_labels) : 24)

    return (GeneratorSet = G,
            BusSet       = B,
            TimeSet      = 1:T,
            PiecewiseSet = 1:j,
            T            = T)
end

end

module Variables
using JuMP
export definir_variables

function definir_variables(modelo::JuMP.Model, set)
    u     = @variable(modelo, u[set.GeneratorSet,set.TimeSet], Bin)
    p     = @variable(modelo, p[set.GeneratorSet,set.TimeSet] >= 0)
    w     = @variable(modelo, w[set.GeneratorSet,set.TimeSet], Bin)
    v     = @variable(modelo, v[set.GeneratorSet,set.TimeSet], Bin)
    theta = @variable(modelo, theta[set.BusSet,set.TimeSet])
    return (u=u, p=p, w=w, v=v, theta=theta)
end

end

module Objetivo
using JuMP
export funcion_objetivo

function funcion_objetivo(modelo::JuMP.Model, par, set, var)
    @objective(modelo, Min,
        sum(par.generators[g].VariableCost  * var.p[g,t] +
             par.generators[g].StartUpCost   * var.u[g,t] +
             par.generators[g].ShutDownCost  * var.v[g,t]
             for g in set.GeneratorSet, t in set.TimeSet)
    )
    return nothing
end

end


