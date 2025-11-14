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



# %% Componente de Cabo

"""Define um componente de cabo coaxial.

Ele é composto por condutor tubular e um isolante externo.

# Atributos

- `radius_in`: Raio interno do condutor [m]. Use zero se sólido.
- `radius_ext`: Raio externo do condutor [m].
- `radius_ext_insulator`: Raio externo do isolante [m].
- `rho_c`: Resistividade elétrica do condutor em Ω/m.
- `mur_c`: Permeabilidade magnética relativa do condutor.
- `mur_d`: Permeabilidade magnética relativa do isolante.
- `epsr`: Permissividade elétrica relativa do isolante.
- `name`: Identificação do objeto.
"""
mutable struct CableComponent
    radius_in::Float64
    radius_ext::Float64
    radius_ext_insulator::Float64
    rho_c::Float64
    mur_c::Float64
    mur_d::Float64
    epsr::Float64
    name::String
    "Atributo privado: índice interno do condutor."
    _index::Union{Nothing, Int}

    # Constructor
    function CableComponent(
        radius_in::Real,
        radius_ext::Real,
        radius_ext_insulator::Real,
        rho_c::Real,
        mur_c::Real,
        mur_d::Real,
        epsr::Real,
        name::String = "componente",
    )
        if radius_in < 0 || radius_ext < 0 || radius_ext_insulator < 0
            throw(ArgumentError("Raio deve ser positivo."))
        end

        if radius_in > radius_ext > radius_ext_insulator
            throw(ArgumentError("Raio externo deve ser maior que o interno."))
        end

        new(
            Float64(radius_in),
            Float64(radius_ext),
            Float64(radius_ext_insulator),
            Float64(rho_c),
            Float64(mur_c),
            Float64(mur_d),
            Float64(epsr),
            String(name),
            nothing
        )
    end
end


function Base.show(io::IO, c::CableComponent)
    println(io, "CableComponent\n",
                "name = $(c.name)\n",
                "radii = [$(c.radius_in*1e3), $(c.radius_ext*1e3), $(c.radius_ext_insulator*1e3)] mm\n",
                "ρ = $(c.rho_c) Ω·m\nμr_c = $(c.mur_c)\nμr_d = $(c.mur_d)\nεr_d = $(c.epsr)"
    )
end


# %% Cabo Coaxial

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


# %% Cabo Pipe-Type (Tubular)


"""Define um cabo Pipe-Type (Tubular).

Ele é composto por isolante interno, condutor, isolante externo,
uma posição (x,y) e outros cables que podem ser ou Pipe ou Coaxial.

# Atributos

- `radius_in`: Raio interno do condutor [m]. Use zero se sólido.
- `radius_ext`: Raio externo do condutor [m].
- `radius_ext_insulator`: Raio externo do isolante [m].
- `rho_c`: Resistividade elétrica do condutor em Ω/m.
- `mur_c`: Permeabilidade magnética relativa do condutor.
- `mur_d_in`: Permeabilidade magnética relativa do isolante interno.
- `mur_d_ext`: Permeabilidade magnética relativa do isolante externo.
- `epsr_in`: Permissividade elétrica relativa do isolante interno.
- `epsr_ext`: Permissividade elétrica relativa do isolante externo.
- `cables`: Lista de cables dentro do Pipe.
- `x`: Posição horizontal [m].
- `y`: Posição vertical [m].
- `name`: Identificação do objeto.
- `sigma_medium`: Condutividade do meio externo [S/m].
- `epsr_medium`: Permissividade elétrica relativa do meio externo.
- `cable_length`: Comprimento do cabo [m].
"""
mutable struct PipeCable <: AbstractCable
    radius_in::Float64
    radius_ext::Float64
    radius_ext_insulator::Float64
    rho_c::Float64
    mur_c::Float64
    mur_d_in::Float64
    mur_d_ext::Float64
    epsr_in::Float64
    epsr_ext::Float64
    cables::Vector{Union{PipeCable, CoaxialCable}}
    x::Float64
    y::Float64
    name::String
    cable_length::Float64
    sigma_medium::Float64
    epsr_medium::Float64
    "Atributo privado: índice interno do condutor."
    _index::Union{Nothing, Int}

    # Constructors
    function PipeCable(
        radius_in::Real,
        radius_ext::Real,
        radius_ext_insulator::Real,
        rho_c::Real,
        mur_c::Real,
        mur_d_in::Real,
        mur_d_ext::Real,
        epsr_in::Real,
        epsr_ext::Real;
        cables = Vector{Union{PipeCable, CoaxialCable}}(),
        x::Real = 0.0,
        y::Real = 0.0,
        name::AbstractString = "Pipe",
        cable_length::Real = 1.0,
        sigma_medium::Real = 5.0,
        epsr_medium::Real = 81.0,
    )
        cobj = deepcopy(collect(cables))

        # Validar raios
        for (i, cabo) in enumerate(cobj)
            r1 = outer_radius(cabo)
            d1 = hypot(cabo.x, cabo.y)
            d2 = d1 + r1
            if d2 > radius_in
                throw(ArgumentError(
                    "O cabo[$(i)], $(cabo.name), está fora do Pipe."
                ))
            end
        end

        new(
            float(radius_in),
            float(radius_ext),
            float(radius_ext_insulator),
            float(rho_c),
            float(mur_c),
            float(mur_d_in),
            float(mur_d_ext),
            float(epsr_in),
            float(epsr_ext),
            cobj,
            float(x),
            float(y),
            String(name),
            float(cable_length),
            float(sigma_medium),
            float(epsr_medium),
            nothing
        )
    end
end


function Base.show(io::IO, pc::PipeCable)
    return println(io, "PipeCable\n",
                       "name = $(pc.name)\n", 
                       "posição = [$(pc.x * 1e3), $(pc.y * 1e3)] mm\n", 
                       "raio_interno = $(pc.radius_in * 1e3) mm\n",
                       "raio_externo = $(pc.radius_ext_insulator * 1e3) mm\n",
                       "cabos_internos = $(length(pc.cables))"
    )
end


function outer_radius(pc::PipeCable)
    return pc.radius_ext_insulator
end


function name_index_dict(pc::PipeCable)
    d = Dict{Union{Int, Missing}, String}()
    if pc._index !== nothing
        d[pc._index] = pc.name
    end
    for c in pc.cables
        for (k, v) in name_index_dict(c)
            d[k] = v
        end
    end
    return d
end


"""Desloca em `(dx, dy)` as coordenadas do cabo e todos seus cabos internos."""
function move_group(pc::PipeCable, dx::Real, dy::Real)
    pc.x += dx
    pc.y += dy
    for c in pc.cables
        if c isa PipeCable
            move_group(c, dx, dy)
        elseif c isa CoaxialCable
            c.x += dx
            c.y += dy
        end
    end
end


# TODO Visualização
