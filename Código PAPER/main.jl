include("modelo/parametros.jl")
include("modelo/variables.jl")
include("modelo/restricciones.jl")

# 2. Incluir el módulo principal que los utiliza
# (Este archivo debe ir al final, ya que usa los de arriba)
include("modelo/unit_commitment.jl") 

# 3. Usar el módulo principal para construir y resolver
using .ModeloUC
using .Restricciones: _bus_key, _F
using DataFrames
using JuMP: value
using MathOptInterface
using XLSX

const MOI = MathOptInterface

function _collect_columns(df::DataFrame)
    names_df = names(df)
    columns = [df[!, name] for name in names_df]
    headers = String.(names_df)
    return columns, headers
end

function export_results_to_excel(par, set, var, filepath::AbstractString)
    mkpath(dirname(filepath))
    rows = NamedTuple{(:generator, :bus_id, :time, :u, :p), Tuple{Int, Any, Int, Float64, Float64}}[]
    for g in set.GeneratorSet, t in set.TimeSet
        push!(rows, (
            generator = g,
            bus_id    = par.generators[g].bus_id,
            time      = t,
            u         = value(var.u[g,t]),
            p         = value(var.p[g,t])
        ))
    end
    gen_df = DataFrame(rows)

    scc_rows = NamedTuple{(:bus_index, :bus_id, :time, :I_syn, :I_ibg, :I_total, :I_limit),
                          Tuple{Int, Any, Int, Float64, Float64, Float64, Float64}}[]
    nb = length(par.buses)
    bus_idx = Dict(_bus_key(par.buses[b].bus_id) => b for b in 1:nb)
    G   = set.GeneratorSet
    C   = hasproperty(set, :IBGSet) ? set.IBGSet : Base.OneTo(0)
    Phi = isempty(C) ? Int[] : [ bus_idx[_bus_key(par.ibgs[c].bus_id)] for c in C ]
    beta = _F(par.scc.beta; default=0.95)
    Vn   = _F(par.scc.Vn; default=1.0)
    Ig = Dict{Int, Float64}()
    for g in G
        xd = _F(par.generators[g].Xdpp; default=0.2)
        Ig[g] = (beta * Vn) / (xd == 0.0 ? 0.2 : xd)
    end
    alpha_val = _F(par.scc.alpha; default=1.0)

    for F in set.BusSet, t in set.TimeSet
        I_syn = 0.0
        for g in G
            I_syn += Ig[g] * value(var.mu[F,g,t])
        end
        I_ibg = 0.0
        if !isempty(C)
            for (c_idx, c) in enumerate(C)
                I_ibg += _F(par.ibgs[c].If_pu; default=1.0) * alpha_val * value(var.Z[F, Phi[c_idx], t])
            end
        end
        I_total = I_syn + I_ibg
        push!(scc_rows, (
            bus_index = F,
            bus_id    = par.buses[F].bus_id,
            time      = t,
            I_syn     = I_syn,
            I_ibg     = I_ibg,
            I_total   = I_total,
            I_limit   = _F(par.buses[F].IminSCC; default=0.0)
        ))
    end
    scc_df = DataFrame(scc_rows)

    gen_columns, gen_headers = _collect_columns(gen_df)
    scc_columns, scc_headers = _collect_columns(scc_df)

    XLSX.writetable(filepath;
        overwrite=true,
        Generadores=(gen_columns, gen_headers),
        Cortocircuito=(scc_columns, scc_headers))
end

println("Construyendo el modelo...")
modelo, par, set, var = ModeloUC.construir_modelo()
bus_map = Dict(_bus_key(par.buses[b].bus_id) => b for b in set.BusSet)
println("Mapa bus_id -> índice interno:")
println(bus_map)
println("Modelo construido. Resolviendo...")
println("Matriz de impedancias:")
impedance_df = DataFrame([(
    from_bus = linea.from_bus,
    to_bus = linea.to_bus,
    R = linea.R,
    X = linea.X,
    B = linea.B
) for linea in par.impedances])
println(impedance_df)

status = ModeloUC.solve_modelo(modelo)

if status == MOI.OPTIMAL
    output_path = joinpath(@__DIR__, "modelo", "data", "output", "resultados.xlsx")
    export_results_to_excel(par, set, var, output_path)
    println("Resultados exportados a: $output_path")
else
    println("No se exportan resultados porque el estado es $status")
end

println("Fin de la ejecución.")
