module Parametros

using DataFrames, XLSX
using Statistics

export read_input_data, GeneratorData, BusData, MatrizImpedancia, IBGData, SCCParams, FreqParams

struct GeneratorData # Generador convencional
    bus_id ::Any
    bus_name ::Any
    type ::Any      # (1: Slack, 2: PV convencional)
    Pmin ::Any
    Pmax ::Any
    VariableCost ::Any
    Ramp ::Any
    StartUpCost ::Any
    ShutDownCost ::Any
    NoLoadCost ::Any
    MinimumUpTime ::Any
    MinimumDownTime ::Any
    H_c::Any
    Xdpp ::Any
end

struct IBGData # Generador en base a inversores
    bus_id ::Any
    bus_name ::Any
    type ::Any      # (3: IBG)
    gamma ::Any
    If_pu  ::Any
    Hs_max ::Any
end

struct BusData
    bus_id ::Any
    bus_name ::Any
    type ::Any      # (0: PQ)
    Pd ::Any
    IminSCC ::Any # Límite fijado (0, 2, 5, 10)
end

struct MatrizImpedancia
    from_bus ::Any
    to_bus ::Any
    R ::Any
    X ::Any
    B ::Any
end

struct SCCParams # Parámetros de cortocircuito
    beta ::Any
    Vn ::Any
    alpha ::Any
end

struct FreqParams
    delta_P ::Any
    Td ::Any
    f_nadir ::Any
    D ::Any
    f_ss ::Any
    f_rocof ::Any
end

_read_df(path::AbstractString, sheet::AbstractString) = DataFrame(XLSX.readtable(path, sheet))

function read_input_data(data_path::String)
    df_bus    = _read_df(data_path, "BUS DATA")
    df_branch = _read_df(data_path, "BRANCH DATA")

    generators = GeneratorData[]
    buses = BusData[]
    impedances = MatrizImpedancia[]
    ibgs = IBGData[]
    scc = SCCParams[]
    freq = FreqParams[]

    for row in eachrow(df_bus)
        demanda = [row.t1, row.t2, row.t3, row.t4, row.t5, row.t6, row.t7, row.t8, row.t9, row.t10,
                   row.t11, row.t12, row.t13, row.t14, row.t15, row.t16, row.t17, row.t18, row.t19, row.t20,
                   row.t21, row.t22, row.t23, row.t24]
        push!(buses, BusData(row.id, row.Nombre, row.Tipo, demanda, 0.0))


        if row.Tipo == 1 || row.Tipo == 2
            push!(generators, GeneratorData(row.id, row.Nombre, row.Tipo,
                                            row.Pmin, row.Pmax, row.C_Variable,
                                            row.Rampa, row.Costo_On, row.Costo_Off,
                                            row.Costo_NoLoad, row.T_min_On, row.T_min_Off, 
                                            row.H, row.X))
        end

        if row.Tipo == 3
            push!(ibgs, IBGData(row.id, row.Nombre, row.Tipo, 0.1, 1, row.H))
        end
    end

    for row in eachrow(df_branch)
        push!(impedances, MatrizImpedancia(row.From, row.To, row.R, row.X, row.B))
    end

    scc = SCCParams(0.95, 1.0, 1)  # alpha = 1.0
    freq = FreqParams(50, 10, 0.8, 0.005, 0.5, 0.5)

    return generators, buses, impedances, ibgs, scc, freq
end

end 
