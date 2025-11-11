module Restricciones
using JuMP
export generar_restricciones

# =========================
#  Funciones Auxiliares
# =========================

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

# Accede a la demanda [t] para el bus [b]
_demand(par, b, t) = _F(par.buses[b].Pd[t]; default=0.0)

# =========================
#  Restricciones Principales
# =========================

function generar_restricciones(modelo, par, set, var; u0=nothing, p0=nothing, enforce_end::Bool=false)
    u, p = var.u, var.p
    y, z = _y(var), _z(var)
    theta = hasproperty(var, :theta) ? var.theta :
            error("Falta `theta[b,t]` en Variables.definir_variables")

    # --- Mapeo de Buses ---
    bus_idx = Dict(_bus_key(par.buses[b].bus_id) => b for b in set.BusSet)
    

    # --- Mapeo de Ramas (DC Flow) ---
    L = length(par.impedances)
    from_idx = Vector{Int}(undef, L)
    to_idx   = Vector{Int}(undef, L)
    Bij      = Vector{Float64}(undef, L)  # 1/X

    for br in 1:L
        f_key = _bus_key(par.impedances[br].from_bus)
        t_key = _bus_key(par.impedances[br].to_bus)

        from_idx[br] = bus_idx[f_key]
        to_idx[br]   = bus_idx[t_key]

        X = _F(par.impedances[br].X; default=0.0)
        Bij[br] = 1.0 / X
    end

    # ========== 1) KCL (balance nodal DC) ==========
    @constraint(modelo, KCL[t in set.TimeSet, b in set.BusSet],
        # Inyección de Generadores Síncronos en el bus b
        sum(p[g,t] for g in set.GeneratorSet if _bus_key(par.generators[g].bus_id) == _bus_key(par.buses[b].bus_id))
        
        # Demanda en el bus b
        - _demand(par, b, t)
        ==
        # Flujo saliente de b
        sum(Bij[br] * (theta[b,t] - theta[to_idx[br],t]) for br in 1:L if from_idx[br] == b) +
        # Flujo entrante a b
        sum(Bij[br] * (theta[b,t] - theta[from_idx[br],t]) for br in 1:L if to_idx[br] == b)
    )

    # ========== 2) Límites Pmin/Pmax (Generadores Síncronos) ==========
    @constraint(modelo, PminLimit[g in set.GeneratorSet, t in set.TimeSet],
        p[g,t] >= _F(par.generators[g].Pmin) * u[g,t]
    )
    @constraint(modelo, PmaxLimit[g in set.GeneratorSet, t in set.TimeSet],
        p[g,t] <= _F(par.generators[g].Pmax) * u[g,t]
    )


    # ========== 3) Lógica on/off y no simultáneo ==========
    if u0 === nothing
        u0_val = 0 # Asumir apagado si no hay estado inicial
        @constraint(modelo, [g in set.GeneratorSet], u[g,1] - u0_val == y[g,1] - z[g,1])
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
        p0_val = 0.0 # Asumir 0 si no hay estado inicial
        @constraint(modelo, [g in set.GeneratorSet],
            p[g,1] - p0_val <= _F(par.generators[g].Ramp) * u[g,1] + _F(par.generators[g].Pmin) * (_y(var))[g,1]
        )
        @constraint(modelo, [g in set.GeneratorSet],
            p0_val - p[g,1] <= _F(par.generators[g].Ramp) * u[g,1] + _F(par.generators[g].Pmin) * (_z(var))[g,1]
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
    
    # --- Restricciones de frecuencia ---
    add_rocof_constraint!(modelo, par, set, var)

    # --- Restricciones de SCC ---
    add_scc_constraints!(modelo, par, set, var)

    return nothing
end

# =========================
#   Frecuencia
# =========================

const DEFAULT_F_NOMINAL = 50.0  # Hz, se usa si no se entrega en los datos

function add_rocof_constraint!(modelo, par, set, var)
    hasproperty(par, :freq) || return

    delta_P = abs(_F(par.freq.delta_P; default=0.0))
    f_rocof = _F(par.freq.f_rocof; default=0.0)
    if delta_P == 0.0 || f_rocof <= 0.0
        return
    end

    f_nominal = if hasproperty(par.freq, :f_base)
        _F(getfield(par.freq, :f_base); default=DEFAULT_F_NOMINAL)
    elseif hasproperty(par.freq, :f_nominal)
        _F(getfield(par.freq, :f_nominal); default=DEFAULT_F_NOMINAL)
    else
        DEFAULT_F_NOMINAL
    end

    inertia_req = (f_nominal * delta_P) / (2.0 * f_rocof)
    inertia_coeff = [ _F(par.generators[g].H_c; default=0.0) * _F(par.generators[g].Pmax; default=0.0)
                      for g in set.GeneratorSet ]
    isempty(inertia_coeff) && return

    u = var.u
    @constraint(modelo, RoCoF[t in set.TimeSet],
        sum(inertia_coeff[g] * u[g,t] for g in set.GeneratorSet) >= inertia_req
    )
end

# =========================
#   SCC
# =========================

function add_scc_constraints!(modelo, par, set, var)
    @assert hasproperty(var, :Z)  "Falta var.Z[i,j,t] en variables.jl"
    @assert hasproperty(var, :mu) "Falta var.mu[i,g,t] en variables.jl"
    @assert hasproperty(var, :u)  "Falta var.u[g,t] (commitment)"

    Z  = var.Z
    mu = var.mu
    u  = var.u

    # --- Mapeos de buses ---
    nb  = length(par.buses)
    bus_idx = Dict(_bus_key(par.buses[b].bus_id) => b for b in 1:nb)

    # --- Ramas y susceptancias (DC: B = 1/X) ---
    L = length(par.impedances)
    from_idx = Vector{Int}(undef, L)
    to_idx   = Vector{Int}(undef, L)
    Bij      = Vector{Float64}(undef, L)

    for ℓ in 1:L
        f_key = _bus_key(par.impedances[ℓ].from_bus)
        t_key = _bus_key(par.impedances[ℓ].to_bus)
        @assert haskey(bus_idx, f_key) "from_bus=$(par.impedances[ℓ].from_bus) no existe en BUS DATA"
        @assert haskey(bus_idx, t_key) "to_bus=$(par.impedances[ℓ].to_bus) no existe en BUS DATA"
        from_idx[ℓ] = bus_idx[f_key]
        to_idx[ℓ]   = bus_idx[t_key]
        X = _F(par.impedances[ℓ].X; default=0.0)
        @assert X != 0.0 "X=0 en la rama $ℓ; el modelo DC requiere X≠0"
        Bij[ℓ] = 1.0 / X
    end

    # --- Y0 (Matriz de admitancia de red pasiva, nb x nb) ---
    Y0 = zeros(Float64, nb, nb)
    for ℓ in 1:L
        i = from_idx[ℓ]; j = to_idx[ℓ]; B = Bij[ℓ]
        Y0[i,i] += B;  Y0[j,j] += B
        Y0[i,j] -= B;  Y0[j,i] -= B
    end

    # --- Índices de generadores síncronos por barra Ψ(g) ---
    G   = set.GeneratorSet
    Psi = [ bus_idx[_bus_key(par.generators[g].bus_id)] for g in G ]

    # --- IBGs (si existen) y sus barras Φ(c) ---
    C   = (hasproperty(set, :IBGSet) ? set.IBGSet : Base.OneTo(0))
    Phi = (isempty(C) ? Int[] : [ bus_idx[_bus_key(par.ibgs[c].bus_id)] for c in C ])

    # --- Corriente Norton SG: I_g = (β·Vn)/Xd″_g ---
    
    # ADAPTACIÓN 1: par.scc.betaE cambiado a par.scc.beta
    beta = _F(par.scc.beta; default=0.95) 
    
    Vn   = _F(par.scc.Vn; default=1.0)
    Xdpp = [ _F(par.generators[g].Xdpp; default=0.2) for g in G ]
    Ig   = [ (beta * Vn) / (Xdpp[g] == 0.0 ? 0.2 : Xdpp[g]) for g in G ] # Evita división por cero

    # --- Acotamos Z y McCormick para μ = Z[:,Ψ(g),t] * u[g,t] ---
    @constraint(modelo, [i=1:nb, j=1:nb, t in set.TimeSet], -_F(par.scc.Z_max; default=10.0) <= Z[i,j,t] <= _F(par.scc.Z_max; default=10.0))
    
    @constraint(modelo, [i=1:nb, g in G, t in set.TimeSet],  mu[i,g,t] <=  Z[i, Psi[g], t] + (1 - u[g,t]) * _F(par.scc.Z_max; default=10.0))
    @constraint(modelo, [i=1:nb, g in G, t in set.TimeSet],  mu[i,g,t] >=  Z[i, Psi[g], t] - (1 - u[g,t]) * _F(par.scc.Z_max; default=10.0))
    @constraint(modelo, [i=1:nb, g in G, t in set.TimeSet],  mu[i,g,t] <=  u[g,t] * _F(par.scc.Z_max; default=10.0))
    @constraint(modelo, [i=1:nb, g in G, t in set.TimeSet],  mu[i,g,t] >= -u[g,t] * _F(par.scc.Z_max; default=10.0))

    # --- Igualdades: Z * (Y0 + Yg(u)) = I ---
    # Yg(u) es una matriz diagonal con 1/Xdpp si el generador g está en el bus j
    @constraint(modelo, [i=1:nb, j=1:nb, t in set.TimeSet],
        sum( Z[i,k,t] * Y0[k,j] for k in 1:nb )
        + sum( (1.0 / Xdpp[g]) * mu[i,g,t] * (j == Psi[g] ? 1.0 : 0.0) for g in G )
        == (i == j ? 1.0 : 0.0)
    )

    # --- Función auxiliar para alpha (parámetro de IBGs) ---
    
    # ADAPTACIÓN 2: Se reemplazó la lógica de Diccionario 
    # para usar el valor 'alpha' simple de tu struct SCCParams.
    function _alpha_c(c::Int, t::Int)
        # Simplemente retorna el valor único de alpha, independientemente de c o t.
        return _F(par.scc.alpha; default=1.0)
    end

    # --- Límite SCC por barra/tiempo ---
    # Ecuación (16a) del paper, reformulada:
    # -sum(I_g * Z_FΨ(g) * u_g) - sum(I_fc * α_c * Z_FΦ(c)) >= I_Flim * Z_FF
    # (Usando mu para Z*u)
    
    @constraint(modelo, [F=1:nb, t in set.TimeSet],
        sum( Ig[g] * mu[F,g,t] for g in G )
        + (isempty(C) ? 0.0 : sum( _F(par.ibgs[c].If_pu; default=1.0) * _alpha_c(c,t) * Z[F, Phi[c], t] for c in C ))
        >= _F(par.buses[F].IminSCC; default=0.0) * Z[F,F,t]
    )

    return nothing
end

end # Fin del Módulo
