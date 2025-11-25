include("modelo/parametros.jl")
include("modelo/variables.jl")
include("modelo/restricciones.jl")
# Nota: No incluimos modelo/unit_commitment.jl porque construiremos el modelo manualmente aquí
# para poder inyectar los parámetros modificados.

using .Parametros
using .Conjuntos
using .Variables
using .Objetivo
using .Restricciones
using .Restricciones: _bus_key, _F

import JuMP
using DataFrames
using JuMP: value, set_optimizer_attribute, compute_conflict!, all_constraints, get_optimizer_attribute, Model
using Gurobi
using MathOptInterface
using XLSX

const MOI = MathOptInterface

# =========================
# FUNCIONES AUXILIARES
# =========================

function _collect_columns(df::DataFrame)
    names_df = names(df)
    columns = [df[!, name] for name in names_df]
    headers = String.(names_df)
    return columns, headers
end

const DEFAULT_FREQUENCY = 50.0

function _freq_nominal(par)
    hasproperty(par, :freq) || return DEFAULT_FREQUENCY
    freq = par.freq
    if hasproperty(freq, :f_base)
        return _F(getfield(freq, :f_base); default=DEFAULT_FREQUENCY)
    elseif hasproperty(freq, :f_nominal)
        return _F(getfield(freq, :f_nominal); default=DEFAULT_FREQUENCY)
    else
        return DEFAULT_FREQUENCY
    end
end

_freq_delta_p(par) = hasproperty(par, :freq) ? abs(_F(par.freq.delta_P; default=0.0)) : 0.0

# =========================
# EXPORTACIÓN DE RESULTADOS
# =========================

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

    # --- REPORTE DE SCC ---
    scc_rows = NamedTuple{(:bus_index, :bus_id, :time, :I_syn, :I_ibg, :I_total, :I_limit_param, :I_limit_real),
                          Tuple{Int, Any, Int, Float64, Float64, Float64, Float64, Float64}}[]
    
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
        term_syn = 0.0
        for g in G
            try term_syn += Ig[g] * value(var.mu[F,g,t]) catch; end
        end
        
        term_ibg = 0.0
        if !isempty(C)
            for (c_idx, c) in enumerate(C)
                try term_ibg += _F(par.ibgs[c].If_pu; default=1.0) * alpha_val * value(var.Z[F, Phi[c_idx], t]) catch; end
            end
        end

        # Calcular valores físicos reales (SCC)
        # IMPORTANTE: Signos negativos para obtener la magnitud física positiva
        I_syn_real   = term_syn
        I_ibg_real   = term_ibg
        I_total_real = (I_syn_real + I_ibg_real) / value(var.Z[F,F,t])  

        # Calcular el límite real (Lado derecho de la restricción)
        I_limit_param = _F(par.buses[F].IminSCC; default=0.0) 
        I_limit_real  = I_limit_param * value(var.Z[F,F,t])  

        push!(scc_rows, (
            bus_index = F,
            bus_id    = par.buses[F].bus_id,
            time      = t,
            I_syn     = I_syn_real,    
            I_ibg     = I_ibg_real,    
            I_total   = I_total_real,  
            I_limit_param = I_limit_param, 
            I_limit_real = I_limit_real   
        ))
    end
    scc_df = DataFrame(scc_rows)

    # --- REPORTE DE FRECUENCIA ---
    freq_sheet = nothing
    if hasproperty(par, :freq)
        freq_rows = NamedTuple{(:time, :inertia_equiv, :rocof, :frequency),
                               Tuple{Int, Float64, Float64, Float64}}[]
        delta_P = _freq_delta_p(par)
        f_nominal = _freq_nominal(par)
        for t in set.TimeSet
            H_sum = 0.0
            for g in set.GeneratorSet
                H_g = _F(par.generators[g].H_c; default=0.0)
                Pmax_g = _F(par.generators[g].Pmax; default=0.0)
                H_sum += H_g * Pmax_g * value(var.u[g,t])
            end
            if H_sum <= 0.0 || delta_P == 0.0
                rocof = 0.0
                freq_val = f_nominal
            else
                rocof = (f_nominal * delta_P) / (2.0 * H_sum)
                freq_val = max(f_nominal - rocof, 0.0)
            end
            push!(freq_rows, (
                time = t,
                inertia_equiv = H_sum,
                rocof = rocof,
                frequency = freq_val
            ))
        end
        freq_df = DataFrame(freq_rows)
        freq_columns, freq_headers = _collect_columns(freq_df)
        freq_sheet = (freq_columns, freq_headers)
    end

    gen_columns, gen_headers = _collect_columns(gen_df)
    scc_columns, scc_headers = _collect_columns(scc_df) 

    if freq_sheet === nothing
        XLSX.writetable(filepath;
            overwrite=true,
            Generadores=(gen_columns, gen_headers),
            Cortocircuito=(scc_columns, scc_headers))
    else
        XLSX.writetable(filepath;
            overwrite=true,
            Generadores=(gen_columns, gen_headers),
            Cortocircuito=(scc_columns, scc_headers),
            Frecuencia=freq_sheet)
    end
end

# ==================================================================================
# FUNCIÓN DE EJECUCIÓN DE UN SOLO ESCENARIO
# ==================================================================================
function run_scenario(input_file::String, scc_req_val::Float64)
    println("\n" * "="^60)
    println(" EJECUTANDO ESCENARIO: $input_file")
    println(" REQUERIMIENTO IminSCC = $scc_req_val p.u.")
    println("="^60)

    # 1. LEER DATOS CRUDOS
    # (No construimos el modelo aún, solo leemos los structs)
    gens, buses_raw, branches, ibgs, scc, freq = Parametros.read_input_data(input_file)
    
    # 2. MODIFICAR EL REQUERIMIENTO (IminSCC) EN LOS BUSES
    # Creamos un nuevo vector de buses con el parámetro actualizado
    buses_mod = Parametros.BusData[]
    
    for bus in buses_raw
        # Creamos una nueva instancia de BusData reemplazando IminSCC
        # Conservamos el resto de atributos, incluyendo el B_shunt que añadiste
        new_bus = Parametros.BusData(
            bus.bus_id,
            bus.bus_name,
            bus.type,
            bus.Pd,
            scc_req_val, # <--- AQUÍ SE APLICA EL CAMBIO
            bus.B_shunt
        )
        push!(buses_mod, new_bus)
    end

    # 3. EMPAQUETAR PARÁMETROS ACTUALIZADOS
    par = (generators = gens, 
           buses      = buses_mod, 
           impedances = branches, 
           ibgs       = ibgs, 
           scc        = scc, 
           freq       = freq)

    # 4. CONSTRUIR MODELO (Manualmente para usar los nuevos params)
    println(">> Construyendo modelo matemático...")
    modelo = Model(Gurobi.Optimizer)
    set = Conjuntos.definir_conjuntos(par)
    var = Variables.definir_variables(modelo, set)
    Objetivo.funcion_objetivo(modelo, par, set, var)
    Restricciones.generar_restricciones(modelo, par, set, var)

    # 5. CONFIGURAR SOLVER
    set_optimizer_attribute(modelo, "DualReductions", 0)
    
    # 6. RESOLVER
    println(">> Resolviendo modelo...")
    status = Restricciones.MOI.get(modelo, Restricciones.MOI.TerminationStatus())
    # Nota: Para resolver usamos optimize! directamente
    JuMP.optimize!(modelo)
    status = JuMP.termination_status(modelo)

    # 7. EXPORTAR
    base_name = replace(basename(input_file), ".xlsx" => "") 
    # Nombre refleja el requerimiento: "ReqSCC"
    output_filename = "resultados_$(base_name)_ReqSCC$(scc_req_val).xlsx"
    output_path = joinpath(@__DIR__, "modelo", "data", "output", output_filename)

    if status == MOI.OPTIMAL
        export_results_to_excel(par, set, var, output_path)
        println(">> [ÉXITO] Costo: $(JuMP.objective_value(modelo))")
        println(">> Resultados en: $output_path")
    elseif status == MOI.INFEASIBLE
        println(">> [FALLO] El modelo es INFACTIBLE.")
        println(">> (Esto es esperado si el requerimiento es muy alto)")
    else
        println(">> [INFO] Estado del solver: $status")
    end
end

# ==================================================================================
# BUCLE PRINCIPAL DE EJECUCIÓN
# ==================================================================================

# Lista de archivos de entrada
input_files = [
    joinpath(@__DIR__, "modelo", "data", "input", "datos_finales30_pen20.xlsx"),
    joinpath(@__DIR__, "modelo", "data", "input", "datos_finales30_pen50.xlsx"),
    joinpath(@__DIR__, "modelo", "data", "input", "datos_finales30_pen90.xlsx")
]

# Lista de valores de REQUERIMIENTO (IminSCC) a probar
# Al subir este valor, el modelo debería volverse más costoso o infactible
scc_limit_values = [0.0, 2.0, 3.0, 5.0]

println("Iniciando Batch Run (Variando Requerimiento SCC)...")
println("Archivos: $(length(input_files))")
println("Niveles de Requerimiento: $(length(scc_limit_values))")

for file in input_files
    if !isfile(file)
        println("\n[ERROR] Archivo no encontrado: $file")
        continue
    end

    for val in scc_limit_values
        run_scenario(file, Float64(val))
    end
end

println("\n" * "="^60)
println(" BATCH RUN COMPLETADO")
println("="^60)
