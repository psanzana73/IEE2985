module Parametros

using DataFrames, XLSX
using Statistics

export read_input_data, GeneratorData, BusData, MatrizImpedancia, dump_all

struct GeneratorData
    bus_id ::Any
    bus_name ::Any
    type ::Any #(1: Slack  2: PV)
    Pmin ::Any
    Pmax ::Any
    VariableCost ::Any
    Ramp ::Any
    StartUpCost ::Any
    ShutDownCost ::Any
    MinimumUpTime ::Any
    MinimumDownTime ::Any
end

struct BusData
    bus_id ::Any
    bus_name ::Any
    type ::Any #(0: PQ)
    Pd ::Any
end

struct MatrizImpedancia
    from_bus ::Any
    to_bus ::Any
    R ::Any
    X ::Any
    B ::Any
end

_read_df(path::AbstractString, sheet::AbstractString) = DataFrame(XLSX.readtable(path, sheet))

function _assert_cols(df::DataFrame, req_syms::Vector{Symbol}, where::AbstractString)
    present = Set(propertynames(df))  # símbolos
    missing = [c for c in req_syms if !(c in present)]
    @assert isempty(missing) "Faltan columnas en '$where': $(missing). Presentes: $(collect(present))"
end

function read_input_data(data_path::String)
    df_bus    = _read_df(data_path, "BUS DATA")
    df_branch = _read_df(data_path, "BRANCH DATA")

    generators = GeneratorData[]
    buses = BusData[]
    impedances = MatrizImpedancia[]

    for row in eachrow(df_bus)
        demanda = [row.t1, row.t2, row.t3, row.t4, row.t5, row.t6, row.t7, row.t8, row.t9, row.t10,
                   row.t11, row.t12, row.t13, row.t14, row.t15, row.t16, row.t17, row.t18, row.t19, row.t20,
                   row.t21, row.t22, row.t23, row.t24]
        push!(buses, BusData(row.id, row.Nombre, row.Tipo, demanda))
        if row.Tipo != 0
            push!(generators, GeneratorData(row.id, row.Nombre, row.Tipo,
                                            row.Pmin, row.Pmax, row.C_Variable,
                                            row.Rampa, row.Costo_On, row.Costo_Off,
                                            row.T_min_On, row.T_min_Off))
        end
    end

    for row in eachrow(df_branch)
        push!(impedances, MatrizImpedancia(row.From, row.To, row.R, row.X, row.B))
    end

    return generators, buses, impedances
end

end
