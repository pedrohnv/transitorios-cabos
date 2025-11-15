# Testa o cálculo dos parâmetros elétricos dos cabos.

using Test
using NPZ

include("../src/linhas_transmissao.jl")

@testset "Cable electric parameters" begin
    # Setup cable: tripolar with sheathed core
    # Material properties
    epsr_c = 3.31     # condutor central
    epsr_b = 2.3      # blindagem, jaqueta da veia
    epsr_a = 10.0     # última camada, envolve a armadura

    rho_c = 1.7241e-8   # Ω·m, condutor central
    rho_b = 2.2e-7      # Ω·m, blindagem
    rho_a = 2.86e-8     # Ω·m, armadura

    mur_c = 1.0     # permeabilidade magnética relativa do condutor central
    mur_b = 1.0     # permeabilidade magnética relativa da blindagem
    mur_a = 300.0   # permeabilidade magnética relativa da armadura

    sig_s = 5.0     # condutividade do mar [S/m]
    eps_s = 81.0     # permissividade relativa do mar

    # Geometria do umbilical
    rc = 1e-3 .* [9.6, 17.054, 18.054, 19.50]   # raios dos cabos coaxias
    ra = 1e-3 .* [48, 59, 65]                   # raios da armadura

    # Posições
    x0, y0 = 0.0, 0.0

    d1 = (rc[end] + 1e-12) / cos(deg2rad(30.0))
    xa = x0
    ya = y0 + d1 * sin(deg2rad(90.0))

    xb = x0 + d1 * cos(deg2rad(210.0))
    yb = y0 + d1 * sin(deg2rad(210.0))

    xc = x0 + d1 * cos(deg2rad(330.0))
    yc = y0 + d1 * sin(deg2rad(330.0))

    # Componentes
    fase_nucleo = CableComponent(0.0, rc[1], rc[2], rho_c, mur_c, 1.0, epsr_c)
    fase_blindagem = CableComponent(rc[2], rc[3], rc[4], rho_b, mur_b, 1.0, epsr_b)

    fase_a = CoaxialCable([fase_nucleo, fase_blindagem], xa, ya)
    fase_b = CoaxialCable([fase_nucleo, fase_blindagem], xb, yb)
    fase_c = CoaxialCable([fase_nucleo, fase_blindagem], xc, yc)

    armadura = PipeCable(
        ra[1],
        ra[2],
        ra[3],
        rho_a,
        mur_a,
        1.0,
        1.0,
        epsr_b,
        epsr_a,
        cables = [fase_a, fase_b, fase_c],
        x = 0.0,
        y = 0.0,
        name = "Pipe",
        sigma_medium = 5.0,
        epsr_medium = 81.0,
        cable_length = 1.0
    )

    # Frequências
    nf = 200
    freq = exp10.(range(0, 9, nf))
    freq_s = freq * 2im * pi


    @testset "Impedância de retorno pelo mar" begin
        z = [cZmar(jw, ra[end], sig_s, eps_s) for jw in freq_s]
        esperado = npzread("test/fixtures/z0_mar.npy")
        @test z ≈ esperado
    end


    @testset "Cabo Coaxial" begin
        z = stack([comp_coaxial_cable_impedance(fase_a, jw) for jw in freq_s])
        esperado = npzread("test/fixtures/coaxial_Z.npy")
        @test z ≈ esperado
    
        p = comp_coaxial_cable_elastance(fase_a)
        esperado = npzread("test/fixtures/coaxial_P.npy")
        @test p ≈ esperado
    end


    @testset "Cabo Pipe-Type" begin
        z = stack([comp_pipe_cable_impedance(armadura, jw) for jw in freq_s])
        esperado = npzread("test/fixtures/pipe_Z.npy")
        @test z ≈ esperado

        p = comp_pipe_cable_elastance(armadura)
        esperado = npzread("test/fixtures/pipe_P.npy")
        @test p ≈ esperado
    end


    @testset "Impedance and Admittance matrices" begin
        zc, yc = zy_cabo(armadura, freq_s, sig_s, eps_s)
        zc_esperado = npzread("test/fixtures/tripolar_Z.npy")
        @test zc ≈ zc_esperado

        yc_esperado = npzread("test/fixtures/tripolar_Y.npy")
        @test yc ≈ yc_esperado
    end

end
