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
