# Funções para linhas de transmissão: Ynodal, modos de propagação, etc.

using LinearAlgebra

include("cabos.jl")
include("cabo_admitancia.jl")
include("cabo_impedancia.jl")


"""Cacula as matrizes Z e Y por unidade de comprimento de um sistema.

# Parâmetros

- `cable_system`: O objeto que representa o sistema de cabos.
- `complex_frequencies`: Frequências angulares complexas no formato `c + jw` \\[rad/s\\].
- `sigma_mar`: Permissividade elétrica do mar \\[S/m\\]. Use `nothing` para ignorar o mar.
- `epsr_mar`: Permissividade elétrica do mar.

# Retorna

- `Z`: Matriz de impedância longitudinal por unidade de comprimento \\[Ω/m\\].
    Dimensões (N condutores, N condutores, K frequências).
- `Y`: Matriz de admitância transversal por unidade de comprimento \\[S/m\\].
    Dimensões (N condutores, N condutores, K frequências).
"""
function zy_cabo(
    cable_system::PipeCable,
    complex_frequencies::AbstractVector{<:Complex{T}},
    sigma_mar::Union{Real, Nothing} = 5.0,
    epsr_mar::Union{Real, Nothing} = 81.0,
) where {T <: Real}
    Y = comp_cable_system_admittance(cable_system, complex_frequencies)
    Z = comp_cable_system_impedance(cable_system, complex_frequencies, sigma_mar, epsr_mar)
    return Z, Y
end


"""Calcula a matriz de admitância nodal `yn` de um cabo.

# Parâmetros

- `Z`: Matriz de impedância longitudinal [Ω/m]. Dimensões (N condutores, N condutores).
- `Y`: Matriz de admitância transversal [S/m]. Dimensões (N condutores, N condutores).
- `comprimento`: Comprimento total [m].

# Retorna

- `yn`: Matriz de admitância nodal [S]. Dimensões (2N condutores, 2N condutores).
"""
function ynodal(
    Z::Matrix{<:Complex{T}}, Y::Matrix{<:Complex{T}}, comprimento::Real
) where {T <: Real}
    # Método 1: Decomposição Modal (Autovalores/Autovetores)
    # evals, evecs = eig(Z * Y)
    # d = sqrt(evals)
    # Tv = evecs
    # d = sqrt(evals)
    # Tvi = inv(Tv)
    # inv_zc_Tv = inv(zc) * Tv
    # hm1 = exp(-d * comprimento)
    # Am1 = d * (1 + hm1^2) / (1 - hm1^2)
    # Bm1 = -2.0 * d * hm1 / (1 - hm1^2)
    # YL1 = inv_zc_Tv * diagm(Am1) * Tvi
    # YL2 = inv_zc_Tv * diagm(Bm1) * Tvi

    # Método 2: Exponencial de Matriz (Função de Transferência)
    ZY = Z * Y
    sqrtzy = sqrt(ZY)
    Yc = inv(Z) * sqrtzy
    H = exp(-comprimento * sqrtzy)
    H2 = H * H
    inv_IH2 = inv(I - H2)
    YL1 = Yc * (I + H2) * inv_IH2
    YL2 = -2.0 * Yc * H * inv_IH2

    yn = [
        YL1 YL2;
        YL2 YL1
    ]
    return yn
end


"""Calcula a matriz de admitância nodal `yn` de um cabo para K frequências.

# Parâmetros

- `Z`: Matriz de impedância longitudinal [Ω/m].
    Dimensões (N condutores, N condutores, K frequências).
- `Y`: Matriz de admitância transversal [S/m].
    Dimensões (N condutores, N condutores, K frequências).
- `comprimento`: Comprimento total [m].

# Retorna

- `yn`: Matriz de admitância nodal [S].
    Dimensões (2N condutores, 2N condutores, K frequências).
"""
function ynodal_array(
    Z::Array{<:Complex{T}, 3}, Y::Array{<:Complex{T}, 3}, comprimento::Real
) where {T <: Real}
    nf = size(Z, 3)
    yn = [ynodal(Z[:, :, f], Y[:, :, f], comprimento) for f in 1:nf]
    yn = stack(yn)  # Empilha ao longo do eixo da frequência
    return yn
end
