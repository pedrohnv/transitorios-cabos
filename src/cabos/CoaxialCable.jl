"""Define um cabo coaxial.

Ele é composto por camadas de ComponenteCabo e uma posição (x,y).

# Atributos

- `components`: Lista de CableComponent que compõe o cabo, do mais interno ao mais externo.
- `x`: Posição horizontal [m].
- `y`: Posição vertical [m].
- `name`: Identificação do objeto.
- `cable_length`: Comprimento do cabo [m].
"""
mutable struct CoaxialCable <: AbstractCable
    components::Vector{CableComponent}
    x::Float64
    y::Float64
    name::String
    cable_length::Float64
    "Atributo privado: índice interno do condutor."
    _index::Union{Nothing, Int}

    # Constructors
    function CoaxialCable(
        components::AbstractVector{<:CableComponent},
        x::Real = 0.0,
        y::Real = 0.0,
        name::AbstractString = "Coaxial",
        cable_length::Real = 1.0
    )
        comps = deepcopy(collect(components))

        # Validar raios
        for i in 2:length(comps)
            r0 = comps[i-1].radius_ext_insulator
            r1 = comps[i].radius_in
            if !isapprox(r0, r1; atol=1e-12)
                throw(ArgumentError("Raio interno do Componente[$i] deve ser igual ao raio externo do isolante do Componente[$(i-1)]."))
            end
        end

        new(
            comps,
            float(x),
            float(y),
            String(name),
            float(cable_length),
            nothing
        )
    end
end


function Base.show(io::IO, c::CoaxialCable)
    println(io, "CoaxialCable\n",
                "name = $(c.name)\n",
                "posição = [$(c.x*1e3), $(c.y*1e3)] mm\n",
                "raio_externo = $(outer_radius(c)*1e3) mm\n",
                "componentes = $(length(c.components))"
    )
end


function outer_radius(c::CoaxialCable)
    return c.components[end].radius_ext_insulator
end


"""Coleta os nomes e _index dos condutores no cabo como um Dicionário."""
function name_index_dict(c::CoaxialCable)
    return Dict{Int, String}(
        comp._index => comp.name for comp in c.components if comp._index !== nothing
    )
end
