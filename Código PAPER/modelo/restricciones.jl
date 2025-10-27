module Restricciones
using JuMP
export generar_restricciones

_F(x; default=0.0) = x === missing ? default :
    x isa Number ? Float64(x) :
    try
        parse(Float64, replace(String(x), "," => "."))
    catch
        default
    end

_I(x; default=0) = x === missing ? default :
    x isa Integer ? Int(x) :
    Int(round(_F(x; default=default)))

function _bus_key(x)::String
    x === missing && return "missing"
    if x isa Integer
        return string(Int(x))
    elseif x isa AbstractFloat
        if isfinite(x) && abs(x - round(x)) < 1e-9
            return string(Int(round(x)))
        else
            return string(x)  # deja el float como string
        end
    else
        s = strip(String(x))
        s = replace(s, "," => ".")
        # intenta parsear número
        v = tryparse(Float64, s)
        if v !== nothing && isfinite(v)
            if abs(v - round(v)) < 1e-9
                return string(Int(round(v)))
            else
                return string(v)
            end
        end
        return lowercase(s)
    end
end

_y(var) = hasproperty(var, :y) ? var.y :
          hasproperty(var, :w) ? var.w :
          error("No encuentro variable de arranque: ni `y` ni `w`.")
_z(var) = hasproperty(var, :z) ? var.z :
          hasproperty(var, :v) ? var.v :
          error("No encuentro variable de parada: ni `z` ni `v`.")
_demand(par, b, t) = _F(par.buses[b].Pd[t]; default=0.0)

function generar_restricciones(modelo, par, set, var; u0=nothing, p0=nothing, enforce_end::Bool=false)
    u, p = var.u, var.p
    y, z = _y(var), _z(var)
    theta = hasproperty(var, :theta) ? var.theta :
            error("Falta `theta[b,t]` en Variables.definir_variables")

    bus_idx = Dict(_bus_key(par.buses[b].bus_id) => b for b in set.BusSet)
    known_bus_keys = collect(keys(bus_idx))

    L = length(par.impedances)
    from_idx = Vector{Int}(undef, L)
    to_idx   = Vector{Int}(undef, L)
    Bij      = Vector{Float64}(undef, L)  # 1/X

    for br in 1:L
        f_key = _bus_key(par.impedances[br].from_bus)
        t_key = _bus_key(par.impedances[br].to_bus)

        @assert haskey(bus_idx, f_key) "from_bus=$(par.impedances[br].from_bus) no está en par.buses. " *
            "Clave normalizada='$f_key'. Conocidos=$(known_bus_keys)"
        @assert haskey(bus_idx, t_key) "to_bus=$(par.impedances[br].to_bus) no está en par.buses. " *
            "Clave normalizada='$t_key'. Conocidos=$(known_bus_keys)"

        from_idx[br] = bus_idx[f_key]
        to_idx[br]   = bus_idx[t_key]

        X = _F(par.impedances[br].X; default=0.0)
        @assert X != 0.0 "X=0 o inválido en rama $(br); DC requiere X≠0"
        Bij[br] = 1.0 / X
    end

    # ========== 1) KCL (balance nodal DC) ==========
    @constraint(modelo, KCL[t in set.TimeSet, b in set.BusSet],
        sum(p[g,t] for g in set.GeneratorSet if _bus_key(par.generators[g].bus_id) == _bus_key(par.buses[b].bus_id))
        - _demand(par, b, t)
        ==
        sum(Bij[br] * (theta[b,t] - theta[to_idx[br],t]) for br in 1:L if from_idx[br] == b) +
        sum(Bij[br] * (theta[b,t] - theta[from_idx[br],t]) for br in 1:L if to_idx[br] == b)
    )

    # ========== 2) Límites Pmin/Pmax ==========
    @constraint(modelo, PminLimit[g in set.GeneratorSet, t in set.TimeSet],
        p[g,t] >= _F(par.generators[g].Pmin) * u[g,t]
    )
    @constraint(modelo, PmaxLimit[g in set.GeneratorSet, t in set.TimeSet],
        p[g,t] <= _F(par.generators[g].Pmax) * u[g,t]
    )


    # ========== 3) Lógica on/off y no simultáneo ==========
    if u0 === nothing
        @constraint(modelo, [g in set.GeneratorSet], u[g,1] - 0 == y[g,1] - z[g,1])
    else
        @constraint(modelo, [g in set.GeneratorSet], u[g,1] - _I(u0[g]) == y[g,1] - z[g,1])
    end
    @constraint(modelo, LogicLink[g in set.GeneratorSet, t in 2:set.T],
        u[g,t] - u[g,t-1] == y[g,t] - z[g,t]
    )
    @constraint(modelo, NoSimultaneousStartStop[g in set.GeneratorSet, t in set.TimeSet],
        y[g,t] + z[g,t] <= 1
    )

    # ========== 4) Mínimos de encendido ==========
    for g in set.GeneratorSet
        UT = _I(par.generators[g].MinimumUpTime; default=0)
        if UT > 0 && set.T - UT + 1 >= 1
            @constraint(modelo, [t in 1:set.T-UT+1],
                sum(u[g,τ] for τ in t:(t+UT-1)) >= UT * y[g,t]
            )
            if enforce_end
                for t in max(1, set.T-UT+2):set.T
                    @constraint(modelo,
                        sum(u[g,τ] for τ in t:set.T) >= (set.T - t + 1) * y[g,t]
                    )
                end
            end
        end
    end

    # ========== 5) Mínimos de apagado ==========
    for g in set.GeneratorSet
        DT = _I(par.generators[g].MinimumDownTime; default=0)
        if DT > 0 && set.T - DT + 1 >= 1
            @constraint(modelo, [t in 1:set.T-DT+1],
                sum(1 - u[g,τ] for τ in t:(t+DT-1)) >= DT * z[g,t]
            )
            if enforce_end
                for t in max(1, set.T-DT+2):set.T
                    @constraint(modelo,
                        sum(1 - u[g,τ] for τ in t:set.T) >= (set.T - t + 1) * z[g,t]
                    )
                end
            end
        end
    end

    # ========== 6) Rampas up/down ==========

    # t = 1
    if p0 === nothing
        @constraint(modelo, [g in set.GeneratorSet],
            p[g,1] - 0.0 <= _F(par.generators[g].Ramp) * u[g,1] + _F(par.generators[g].Pmin) * (_y(var))[g,1]
        )
        @constraint(modelo, [g in set.GeneratorSet],
            0.0 - p[g,1] <= _F(par.generators[g].Ramp) * u[g,1] + _F(par.generators[g].Pmin) * (_z(var))[g,1]
        )
    else
        @constraint(modelo, [g in set.GeneratorSet],
            p[g,1] - _F(p0[g]; default=0.0) <= _F(par.generators[g].Ramp) * u[g,1] + _F(par.generators[g].Pmin) * (_y(var))[g,1]
        )
        @constraint(modelo, [g in set.GeneratorSet],
            _F(p0[g]; default=0.0) - p[g,1] <= _F(par.generators[g].Ramp) * u[g,1] + _F(par.generators[g].Pmin) * (_z(var))[g,1]
        )
    end

    # t >= 2
    @constraint(modelo, [g in set.GeneratorSet, t in 2:set.T],
        p[g,t] - p[g,t-1] <= _F(par.generators[g].Ramp) * u[g,t-1] + _F(par.generators[g].Pmin) * (_y(var))[g,t]
    )
    @constraint(modelo, [g in set.GeneratorSet, t in 2:set.T],
        p[g,t-1] - p[g,t] <= _F(par.generators[g].Ramp) * u[g,t] + _F(par.generators[g].Pmin) * (_z(var))[g,t]
    )

    return nothing

    end
end
