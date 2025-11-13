# Estruturas que definem cabos

using PlotlyJS


# %% Componente de Cabo

"""Define um componente de cabo coaxial.

Ele é composto por condutor tubular e um isolante externo.

Atributos
---------
radius_in: float
    Raio interno do condutor [m]. Use zero se sólido.
radius_ext: float
    Raio externo do condutor [m].
radius_ext_insulator: float
    Raio externo do isolante [m].
rho_c: float
    Resistividade elétrica do condutor em Ω/m.
mur_c: float
    Permeabilidade magnética relativa do condutor.
mur_d: float
    Permeabilidade magnética relativa do isolante.
epsr: float
    Permissividade elétrica relativa do isolante.
name: str
    Identificação do objeto.
"""
mutable struct CableComponent
    "Raio interno do condutor [m]. Use zero se sólido."
    radius_in::Float64

    "Raio externo do condutor [m]."
    radius_ext::Float64

    "Raio externo do isolante [m]."
    radius_ext_insulator::Float64

    "Resistividade elétrica do condutor em Ω/m."
    rho_c::Float64

    "Permeabilidade magnética relativa do condutor."
    mur_c::Float64

    "Permeabilidade magnética relativa do isolante."
    mur_d::Float64

    "Permissividade elétrica relativa do isolante."
    epsr::Float64

    "Identificação do objeto."
    name::String

    "Atributo privado: índice interno do condutor."
    _index::Union{Nothing, Int}

    # Constructors
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


"""Dicionário do CableCombonent."""
function to_dict(c::CableComponent)
    d = Dict{String, Any}()
    for k in fieldnames(CableComponent)
        d[string(k)] = getfield(c, k)
    end
    d["id"] = "CableComponent"
    return deepcopy(d)
end


"""Reconstrói um CableComponent a partir de um Dicionário."""
function cablecomponent_from_dict(d::Dict{String, Any})
    if get(d, "id", "") != "CableComponent"
        throw(ArgumentError("dicionário não representa um CableComponent"))
    end

    return CableComponent(
        d["radius_in"],
        d["radius_ext"],
        d["radius_ext_insulator"],
        d["rho_c"],
        d["mur_c"],
        d["mur_d"],
        d["epsr"];
        name = d["name"]
    )
end


# %% Cabo Coaxial

"""Define um cabo coaxial.

Ele é composto por camadas de ComponenteCabo e uma posição (x,y).

Atributos
---------
components: ComponenteCabo
    Lista de components que compõe o cabo, do mais interno ao externo.
x: float
    Posição horizontal [m].
y: float
    Posição vertical [m].
name: str
    Identificação do objeto.
sigma_medium: float
    Condutividade do meio externo [S/m].
epsr_medium: float
    Permissividade elétrica relativa do meio externo.
comprimento: float
    Comprimento do cabo [m].
"""
mutable struct CoaxialCable
    "Lista de components que compõe o cabo, do mais interno ao externo."
    components::Vector{CableComponent}

    "Posição horizontal [m]."
    x::Float64

    "Posição vertical [m]."
    y::Float64

    "Identificação do objeto."
    name::String

    "Condutividade do meio externo [S/m]."
    sigma_medium::Float64

    "Permissividade elétrica relativa do meio externo."
    epsr_medium::Float64

    "Comprimento do cabo [m]."
    comprimento::Float64

    "Atributo privado: índice interno do condutor."
    _index::Union{Nothing, Int}

end


# Constructors
function CoaxialCable(
    components::AbstractVector{<:CableComponent},
    x::Real = 0.0,
    y::Real = 0.0,
    name::AbstractString = "Coaxial",
    sigma_medium::Real = 81.0,
    epsr_medium::Real = 5.0,
    comprimento::Real = 1.0
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

    CoaxialCable(
        comps,
        float(x),
        float(y),
        String(name),
        float(sigma_medium),
        float(epsr_medium),
        float(comprimento),
        nothing
    )
end


function Base.show(io::IO, c::CoaxialCable)
    println(io, "CoaxialCable\n",
                "name = $(c.name)\n",
                "posição = [$(c.x*1e3), $(c.y*1e3)] mm\n",
                "raio_externo = $(outer_radius(c)*1e3) mm\n",
                "componentes = $(length(c.components))"
    )
end


"""Retorna o raio externo do cabo."""
function outer_radius(c::CoaxialCable)
    return c.components[end].radius_ext_insulator
end


"""Dicionário do CoaxialCable."""
function to_dict(c::CoaxialCable)
    d = Dict{String, Any}()
    for k in fieldnames(CoaxialCable)
        d[string(k)] = getfield(c, k)
    end
    d["id"] = "CoaxialCable"
    return deepcopy(d)
end


"""Reconstrói um CoaxialCable a partir de um Dicionário."""
function coaxialcable_from_dict(d::Dict{String, Any})
    if get(d, "id", nothing) != "CoaxialCable"
        throw(ArgumentError("dicionário não representa um CoaxialCable"))
    end

    comp = [CableComponent.from_dict(cd) for cd in d["components"]]
    CoaxialCable(
        comp;
        x = d["x"],
        y = d["y"],
        name = d["name"],
        sigma_medium = d["sigma_medium"],
        epsr_medium = d["epsr_medium"],
        comprimento = d["comprimento"]
    )
end


"""Coleta os nomes e _index dos condutores no CoaxialCable como um Dicionário."""
function name_index_dict(c::CoaxialCable)
    # Assumes _index on components may be Nothing or Int; skip Nones
    return Dict{Int, String}(
        comp._index => comp.name for comp in c.components if comp._index !== nothing
    )
end


# %% Cabo Tubular Oco


"""Define um cabo Pipe-Type.

Ele é composto por isolante interno, condutor, isolante externo,
uma posição (x,y) e outros cables que podem ser ou Pipe ou Coaxial.

Atributos
---------
radius_in: float
    Raio interno do condutor [m]. Use zero se sólido.
radius_ext: float
    Raio externo do condutor [m].
radius_ext_insulator: float
    Raio externo do isolante [m].
rho_c: float
    Resistividade elétrica do condutor em Ω/m.
mur_c: float
    Permeabilidade magnética relativa do condutor.
mur_d_in: float
    Permeabilidade magnética relativa do isolante interno.
mur_d_ext: float
    Permeabilidade magnética relativa do isolante externo.
epsr_in: float
    Permissividade elétrica relativa do isolante interno.
epsr_ext: float
    Permissividade elétrica relativa do isolante externo.
cables:
    Lista de cables dentro do Pipe.
x: float
    Posição horizontal [m].
y: float
    Posição vertical [m].
name: str
    Identificação do objeto.
sigma_medium: float
    Condutividade do meio externo [S/m].
epsr_medium: float
    Permissividade elétrica relativa do meio externo.
comprimento: float
    Comprimento do cabo [m].
"""
mutable struct PipeCable
    "Raio interno do condutor [m]. Use zero se sólido."
    radius_in::Float64

    "Raio externo do condutor [m]."
    radius_ext::Float64

    "Raio externo do isolante [m]."
    radius_ext_insulator::Float64

    "Resistividade elétrica do condutor em Ω/m."
    rho_c::Float64

    "Permeabilidade magnética relativa do condutor."
    mur_c::Float64

    "Permeabilidade magnética relativa do isolante interno."
    mur_d_in::Float64

    "Permeabilidade magnética relativa do isolante externo."
    mur_d_ext::Float64

    "Permissividade elétrica relativa do isolante interno."
    epsr_in::Float64

    "Permissividade elétrica relativa do isolante externo."
    epsr_ext::Float64

    "Lista de cables dentro do Pipe."
    cables::Vector{Union{PipeCable, CoaxialCable}}

    "Posição horizontal [m]."
    x::Float64

    "Posição vertical [m]."
    y::Float64

    "Identificação do objeto."
    name::String

    "Condutividade do meio externo [S/m]."
    sigma_medium::Float64

    "Permissividade elétrica relativa do meio externo."
    epsr_medium::Float64

    "Comprimento do cabo [m]."
    comprimento::Float64

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
        sigma_medium::Real = 5.0,
        epsr_medium::Real = 81.0,
        comprimento::Real = 1.0
    )
        cobj = deepcopy(collect(cables))

        # Validar raios
        for (i, cabo) in pairs(cobj)
            if i > 1
                r1 = outer_radius(cabo)
                d1 = hypot(cabo.x, cabo.y)
                d2 = d1 + r1
                if d2 > radius_in
                    throw(ValueError(
                        "O cabo[$(i)], $(cabo.name), está fora do Pipe."
                    ))
                end
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
            float(sigma_medium),
            float(epsr_medium),
            float(comprimento),
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


function to_dict(pc::PipeCable)
    d = Dict{String, Any}(
        "radius_in" => pc.radius_in,
        "radius_ext" => pc.radius_ext,
        "radius_ext_insulator" => pc.radius_ext_insulator,
        "rho_c" => pc.rho_c,
        "mur_c" => pc.mur_c,
        "mur_d_in" => pc.mur_d_in,
        "mur_d_ext" => pc.mur_d_ext,
        "epsr_in" => pc.epsr_in,
        "epsr_ext" => pc.epsr_ext,
        "x" => pc.x,
        "y" => pc.y,
        "name" => pc.name,
        "sigma_medium" => pc.sigma_medium,
        "epsr_medium" => pc.epsr_medium,
        "comprimento" => pc.comprimento,
        "_index" => pc._index,
        "cables" => [to_dict(c) for c in pc.cables]
    )
    d["id"] = "PipeCable"
    return d
end


"""Reconstrói um PipeCable a partir de um Dicionário."""
function pipecable_from_dict(d::Dict{String, Any})
    if get(d, "id", nothing) != "PipeCable"
        throw(ArgumentError("dicionário não representa um PipeCable"))
    end

    cables = [c == nothing ? nothing : (c["id"] == "CoaxialCable" ? CoaxialCable.from_dict(c) :
                                      c["id"] == "PipeCable" ? PipeCable.from_dict(c) :
                                      error("tipo de cabo não reconhecido: $(c["id"])")) for c in d["cables"]]

    # Filter out possible nothings and construct PipeCable
    pipes = [cc for cc in cables if cc !== nothing]

    PipeCable(
        d["radius_in"],
        d["radius_ext"],
        d["radius_ext_insulator"],
        d["rho_c"],
        d["mur_c"],
        d["mur_d_in"],
        d["mur_d_ext"],
        d["epsr_in"],
        d["epsr_ext"];
        cables = pipes,
        x = d["x"],
        y = d["y"],
        name = d["name"],
        sigma_medium = d["sigma_medium"],
        epsr_medium = d["epsr_medium"],
        comprimento = d["comprimento"]
    )
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
    pc.x += float(dx)
    pc.y += float(dy)
    for c in pc.cables
        if typeof(c) <: PipeCable
            move_group(c, dx, dy)
        elseif typeof(c) <: CoaxialCable
            c.x += float(dx)
            c.y += float(dy)
        end
    end
end


"""Conta o número de condutores dentro de um cabo."""
function count_conductors_cable(cable)
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


# TODO Visualização
