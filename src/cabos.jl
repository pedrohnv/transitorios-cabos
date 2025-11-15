# Estruturas que definem cabos

"""Dicionário do objeto."""
function struct_to_dict(s)
    d = Dict{String, Any}()
    for key in propertynames(s)
        # Skip private fields (those starting with _)
        if !startswith(string(key), "_")
            d[string(key)] = getfield(s, key)
        end
    end
    d["tipo"] = string(typeof(s))
    return d
end


"""Reconstrói um objeto a partir de um Dicionário."""
function struct_from_dict(d::Dict{String, Any})
    tipo = d["tipo"]
    T = eval(Symbol(tipo))
    #vals = [d[string(f)] for f in fieldnames(T)]
    vals = []
    for f in fieldnames(T)
        # Skip private fields (those starting with _)
        key = string(f)
        if !startswith(key, "_")
            push!(vals, d[key])
        end
    end
    c = T(vals...)
    return c
end


"""Um cabo abstrato deve conter.

# Atributos

- `x`: Posição horizontal [m].
- `y`: Posição vertical [m].
- `name`: Identificação do objeto.
- `cable_length`: Comprimento do cabo [m].
- `_index`: Atributo privado: índice interno do condutor.
"""
abstract type AbstractCable end


"""Retorna o raio externo do cabo."""
function outer_radius(c::AbstractCable)
    return nothing
end


"""Conta o número de condutores dentro de um cabo."""
function count_conductors_cable(cable::AbstractCable)
    if cable isa CoaxialCable
        return length(cable.components)
    elseif cable isa PipeCable
        count = 1
        for inner in cable.cables
            count += count_conductors_cable(inner)
        end
        return count
    else
        return 0
    end
end


include("cabos/CableComponent.jl")
include("cabos/CoaxialCable.jl")
include("cabos/PipeCable.jl")


# TODO Visualização
