#=
Este arquivo executa os testes assertivos dos códigos.
=#

include("param_cabo.jl")

import Combinatorics.combinations

# %% TODO Cálculo dos parâmetros elétricos

function teste_modelo_cabo()
    L = 2.5e3  # m
    # Geometria e parâmetros do umbilical
    # Geometria do umbilical PT
    rc = 1e-3 * ([9.6, 17.054, 18.054, 19.50])  # raios dos cabos SC [m]: condutor central, isolante, blindagem e isolante externo
    ra = 1e-3 * ([48, 59, 65])  # raios da armadura [m]: interno, externo, capa
    nc = 7  # número de condutores (fase + blindagem + armadura)
    nx = 2  # número de condutores "ativos" em cada cabo SC (single core)

    # distância do centro do umbilical ao centro das veias
    d01 = rc[4] / cos(deg2rad(30))
    d02 = d01
    di = [d01, d02, d02]
    theta_jk = (2 * pi) / 3  # ângulo entre cabos do umbilical

    # permissividades elétricas
    epsr_c = 3.31  # condutor central
    epsr_b = 2.3  # blindagem, jaqueta da veia
    epsr_a = 10  # última camada, que envolve a armadura

    # resistividades elétricas [Ohm.m]
    rho_c = 1.7241e-8  # condutor central
    rho_b = 2.2e-7  # blindagem, jaqueta da veia
    rho_a = 2.86e-8  # armadura

    mur_a = 300  # permeabilidade magnética relativa da armadura
    sig_s = 5.0  # condutividade elétrica do mar [S/m]
    eps_s = 81.0  # permissividade elétrica relativa do mar
    
    nf = 200
    freq = exp10.(range(0, 9, length=nf))
    
    Yc = zeros(ComplexF64, (nf, nc, nc))
    Tv = zeros(ComplexF64, (nf, nc, nc))
    gama = zeros(ComplexF64, (nf, nc))
    hm = zeros(ComplexF64, (nf, nc))
    evelho = zeros(ComplexF64, (nc, nc))
    evals = zeros(ComplexF64, nc)
    for linha = 1:nf
        w = 2 * pi * freq[linha]
        Zcabo, Ycabo = zy_cabo(w, rc, ra, rho_a, rho_b, rho_c,
                               epsr_a, epsr_b, epsr_c, mur_a,
                               theta_jk, di, sig_s, eps_s, nx)
        ZY = Zcabo * Ycabo
        YZ = Ycabo * Zcabo
        Yc[linha, :, :] = inv(Zcabo) * sqrt(ZY)
        if linha == 1
            evals[:], evecs = eigen(YZ)  # right eigenvecs
            # evecs transposto em relação ao Mathematica
            evelho[:] = rot(evecs)  # linha modificada
        else
            evals[:], evelho[:] = eigNR(YZ, evelho, evals, 1e-5, 20000)
        end
        Tv[linha,:,:] = evelho
        d = sqrt.(evals)
        hm[linha,:] = exp.(-d * L)
        gama[linha,:] = d
    end
    p = plot(freq, real.(gama), xaxis=:log, yaxis=:log, label="",
             xlabel="frequência [Hz]", ylabel="Re(γ)")
    display(p)
    # ----------------------
    v = zeros((nf, nc))
    for i = 1:nc
        v[:, i] .= @. 2e-6 * pi * freq / imag(gama[:, i])
    end
    p = plot(freq, v, xaxis=:log, yaxis=:log, label="",
             xlabel="frequência [Hz]", ylabel="v [m/μs]")
    display(p)
    # ----------------------
    p = plot(freq, abs.(hm), xaxis=:log, label="",
             xlabel="frequência [Hz]", ylabel="Amplitude [p.u.]")
    display(p)
end
#using Plots
#teste_modelo_cabo()


# %% Montar YN

function testar_montar_yn_sem_falha()
    # Parâmetros do cabo
    fases = [1]
    blindagens = [2]
    armadura = 3
    condutores = [fases; blindagens; armadura]
    nc1 = length(condutores)
    nc2 = 2 * nc1
    pares_falhas = [
        [1, 2],
        [1, 3],
        [2, 3]
    ]
    num_pares_falhas = length(pares_falhas)  # Int(nc1 * (nc1 - 1) // 2)
    num_segmentos = 1
    segmento_falha = 0
    terminal_fonte = 1
    aterrar_receptor = [true]

    # YN correto (esperado)
    nsis = nc1 * (num_segmentos + 1)
    yn_correto = zeros(Int, nsis, nsis)
    
    yn_correto[terminal_fonte, terminal_fonte] += 1

    # y_emissor_shunt:
    for i in 1:nc1
        yn_correto[i,i] += 1
    end

    # primeiro segmento:
    m1 = 1:(2*nc1)
    for i in m1
        for k in m1
            yn_correto[i,k] += 1
        end
    end
    
    # receptor:
    k2 = nc1 + armadura
    yn_correto[k2, k2] += 1  # armadura para terra
    # blindagens para armadura
    for i in blindagens
        k1 = nc1 + i
        yn_correto[k1, k1] += 1
        yn_correto[k1, k2] += 1
        yn_correto[k2, k1] += 1
        yn_correto[k2, k2] += 1
    end

    yn_correto2 = [
        3 1 1 1 1 1 
        1 2 1 1 1 1 
        1 1 2 1 1 1 
        1 1 1 1 1 1 
        1 1 1 1 2 2 
        1 1 1 1 2 3 
    ]
    @assert yn_correto == yn_correto2

    # montar YN com valores unitários
    yn_segmentos = -ones(ComplexF64, (nc2, nc2, num_segmentos))
    for i = 1:num_segmentos
        for k = 1:nc2
            yn_segmentos[k,k,i] *= -1
        end
    end

    y_curto = 1.0 + 0.0im
    y_fonte = 1.0 + 0.0im
    y_emissor_shunt = ones(ComplexF64, nc1)
    y_falha_serie = zeros(ComplexF64, nc1)
    y_falha_shunt = zeros(ComplexF64, num_pares_falhas)
    yn = montar_yn(yn_segmentos, y_curto, y_fonte, y_emissor_shunt,
                   y_falha_shunt, y_falha_serie, segmento_falha, pares_falhas,
                   terminal_fonte, fases, blindagens, armadura, aterrar_receptor)
    yn = Int.(abs.(yn))
    return all(isapprox.(yn - yn_correto, 0))
end
@assert testar_montar_yn_sem_falha()


function testar_montar_yn_falha_shunt()
    # Parâmetros do cabo
    fases = [1]
    blindagens = [2]
    armadura = 3
    condutores = [fases; blindagens; armadura]
    nc1 = length(condutores)
    nc2 = 2 * nc1
    pares_falhas = [
        [1, 2],
        [1, 3],
        [2, 3]
    ]
    num_pares_falhas = length(pares_falhas)  # Int(nc1 * (nc1 - 1) // 2)
    num_segmentos = 2
    segmento_falha = 1
    terminal_fonte = 1
    aterrar_receptor = [false, true]

    # YN correto (esperado)
    nsis = nc1 * (num_segmentos + 1)
    no_fonte = nsis - 1
    yn_correto = zeros(Int, nsis, nsis)
    
    yn_correto[terminal_fonte, terminal_fonte] += 1
    
    # y_emissor_shunt:
    for i in 1:nc1
        yn_correto[i,i] += 1
    end

    # primeiro segmento:
    m1 = 1:(2*nc1)
    for i in m1
        for k in m1
            yn_correto[i,k] += 1
        end
    end

    # segundo segmento:
    m1 = (nc1+1):(3*nc1)
    for i in m1
        for k in m1
            yn_correto[i,k] += 1
        end
    end

    # y_falhas_shunt:
    m1 = (nc1+1):(2*nc1)
    for (i,k) in pares_falhas
        i += nc1
        k += nc1
        yn_correto[i,i] += 1
        yn_correto[i,k] += 1
        yn_correto[k,i] += 1
        yn_correto[k,k] += 1
    end

    # receptor:
    k2 = 2*nc1 + armadura
    yn_correto[k2, k2] += 1  # armadura para terra
    # blindagens para armadura
    for i in blindagens
        k1 = 2*nc1 + i
        yn_correto[k1, k1] += 1
        yn_correto[k1, k2] += 1
        yn_correto[k2, k1] += 1
        yn_correto[k2, k2] += 1
    end

    yn_correto2 = [
        3 1 1 1 1 1 0 0 0 
        1 2 1 1 1 1 0 0 0 
        1 1 2 1 1 1 0 0 0 
        1 1 1 4 3 3 1 1 1 
        1 1 1 3 4 3 1 1 1 
        1 1 1 3 3 4 1 1 1 
        0 0 0 1 1 1 1 1 1 
        0 0 0 1 1 1 1 2 2 
        0 0 0 1 1 1 1 2 3  
    ]
    @assert yn_correto == yn_correto2

    # montar YN com valores unitários
    yn_segmentos = -ones(ComplexF64, (nc2, nc2, num_segmentos))
    for i = 1:num_segmentos
        for k = 1:nc2
            yn_segmentos[k,k,i] *= -1
        end
    end

    y_curto = 1.0 + 0.0im
    y_fonte = 1.0 + 0.0im
    y_emissor_shunt = ones(ComplexF64, nc1)
    y_falha_serie = zeros(ComplexF64, nc1)
    y_falha_shunt = ones(ComplexF64, num_pares_falhas)
    yn = montar_yn(yn_segmentos, y_curto, y_fonte, y_emissor_shunt,
                   y_falha_shunt, y_falha_serie, segmento_falha, pares_falhas,
                   terminal_fonte, fases, blindagens, armadura, aterrar_receptor)
    yn = Int.(abs.(yn))
    return all(isapprox.(yn - yn_correto, 0))
end
@assert testar_montar_yn_falha_shunt()


function testar_montar_yn_falha_serie()
    # Parâmetros do cabo
    fases = [1]
    blindagens = [2]
    armadura = 3
    condutores = [fases; blindagens; armadura]
    nc1 = length(condutores)
    nc2 = 2 * nc1
    pares_falhas = [
        [1, 2],
        [1, 3],
        [2, 3]
    ]
    num_pares_falhas = length(pares_falhas)  # Int(nc1 * (nc1 - 1) // 2)
    num_segmentos = 2
    segmento_falha = 1
    terminal_fonte = 1
    aterrar_receptor = [false, true]

    # YN correto (esperado)
    nsis = nc1 * (num_segmentos + 2)
    yn_correto = zeros(Int, nsis, nsis)
    
    yn_correto[terminal_fonte, terminal_fonte] += 1
    
    # y_emissor_shunt:
    for i in 1:nc1
        yn_correto[i,i] += 1
    end

    # primeiro segmento:
    m1 = 1:(2*nc1)
    for i in m1
        for k in m1
            yn_correto[i,k] += 1
        end
    end

    # y_falhas_shunt:
    m1 = (nc1+1):(2*nc1)
    for (i,k) in pares_falhas
        i += nc1
        k += nc1
        yn_correto[i,i] += 1
        yn_correto[i,k] += 1
        yn_correto[k,i] += 1
        yn_correto[k,k] += 1
    end

    # y_falhas_serie:
    m1 = (nc1+1):(2*nc1)
    for i = 1:nc1
        i1 = m1[i]
        i2 = m1[i] + nc1
        yn_correto[i1, i1] += 1
        yn_correto[i1, i2] += 1
        yn_correto[i2, i1] += 1
        yn_correto[i2, i2] += 1
    end

    # segundo segmento:
    m1 = (2*nc1+1):(4*nc1)
    for i in m1
        for k in m1
            yn_correto[i,k] += 1
        end
    end

    # receptor:
    k2 = 3*nc1 + armadura
    yn_correto[k2, k2] += 1  # armadura para terra
    # blindagens para armadura
    for i in blindagens
        k1 = 3*nc1 + i
        yn_correto[k1, k1] += 1
        yn_correto[k1, k2] += 1
        yn_correto[k2, k1] += 1
        yn_correto[k2, k2] += 1
    end
    
    yn_correto2 = [
        3 1 1 1 1 1 0 0 0 0 0 0 
        1 2 1 1 1 1 0 0 0 0 0 0 
        1 1 2 1 1 1 0 0 0 0 0 0 
        1 1 1 4 2 2 1 0 0 0 0 0 
        1 1 1 2 4 2 0 1 0 0 0 0 
        1 1 1 2 2 4 0 0 1 0 0 0 
        0 0 0 1 0 0 2 1 1 1 1 1 
        0 0 0 0 1 0 1 2 1 1 1 1 
        0 0 0 0 0 1 1 1 2 1 1 1 
        0 0 0 0 0 0 1 1 1 1 1 1 
        0 0 0 0 0 0 1 1 1 1 2 2 
        0 0 0 0 0 0 1 1 1 1 2 3  
    ]
    @assert yn_correto == yn_correto2

    # montar YN com valores unitários
    yn_segmentos = -ones(ComplexF64, (nc2, nc2, num_segmentos))
    for i = 1:num_segmentos
        for k = 1:nc2
            yn_segmentos[k,k,i] *= -1
        end
    end

    y_curto = 1.0 + 0.0im
    y_fonte = 1.0 + 0.0im
    y_emissor_shunt = ones(ComplexF64, nc1)
    y_falha_serie = ones(ComplexF64, nc1)
    y_falha_shunt = ones(ComplexF64, num_pares_falhas)
    yn = montar_yn(yn_segmentos, y_curto, y_fonte, y_emissor_shunt,
                   y_falha_shunt, y_falha_serie, segmento_falha, pares_falhas,
                   terminal_fonte, fases, blindagens, armadura, aterrar_receptor)
    yn = Int.(abs.(yn))
    return all(isapprox.(yn - yn_correto, 0))
end
@assert testar_montar_yn_falha_serie()

println("Testes para montar_ynodal bem sucedidos.")


# %% Simulação do cabo

begin  # Parâmetros
    ## Tempo
    dt = 2.5e-9  # passo no tempo [s]
    # Usar número nt que seja potência de 2, para diminuir lixo numérico do zero-padding da FFT
    nt = 1024 * 98  # número de passos no tempo
    tempo = range(0, step=dt, length=nt)
    tmax = tempo[end]
    tempo_abre = 10e-6  # [s] que a chave desconecta a fonte do cabo
    nt_trunc = Int(ceil(nt * 0.80))  # 20% dos tempos finais é lixo numérico
    tmax_trunc = tempo[nt_trunc]  # tempo [s] final que será salvo

    ## Geometria do umbilical PT
    Lc = 2.5e3  # comprimento total [m]
    rc = 1e-3 * ([9.6, 17.054, 18.054, 19.50])  # raios dos cabos SC [m]: condutor central, isolante, blindagem e isolante externo
    ra = 1e-3 * ([48, 59, 65])  # raios da armadura [m]: interno, externo, capa
    nx = 2  # número de condutores "ativos" em cada cabo SC (single core)

    # distância do centro do umbilical ao centro das veias
    d01 = rc[4] / cos(deg2rad(30))
    d02 = d01
    di = [d01, d02, d02]
    theta_jk = (2 * pi) / 3  # ângulo entre cabos do umbilical

    # permissividades elétricas
    epsr_c = 3.31  # condutor central
    epsr_b = 2.3  # blindagem, jaqueta da veia
    epsr_a = 10  # última camada, que envolve a armadura

    # resistividades elétricas [Ohm.m]
    rho_c = 1.7241e-8  # condutor central
    rho_b = 2.2e-7  # blindagem, jaqueta da veia
    rho_a = 2.86e-8  # armadura

    mur_a = 300  # permeabilidade magnética relativa da armadura
    sig_s = 5.0  # condutividade elétrica do mar [S/m]
    eps_s = 81.0  # permissividade elétrica relativa do mar

    fases = [1, 3, 5]
    blindagens = [2, 4, 6]
    armadura = 7
    condutores = [fases; blindagens; armadura]
    nc1 = length(condutores)  # número de condutores
    nc2 = nc1 * 2
    nome_condutores = sort(Dict([
        [fases[i] => "fase_$(i)" for i in eachindex(fases)];
        [blindagens[i] => "blindagem_$(i)" for i in eachindex(blindagens)];
        [armadura => "armadura"]
    ]))
    pares_falhas = collect(combinations(condutores, 2))
    num_pares_falhas = length(pares_falhas)

    ## Cálculo dos parâmetros e matrizes
    comprimentos = [0.5, 1.0] .* Lc
    arquivo_matrizes = "teste/umbilical_tripolar.h5"
    precalc_matrizes(arquivo_matrizes, comprimentos, tmax, nt,
                     rc, ra, rho_a, rho_b, rho_c, epsr_a, epsr_b, epsr_c,
                     mur_a, theta_jk, di, sig_s, eps_s, nx)
    #
    v_fonte_t = ones(nt)
    R_fonte = 50.0
    R_curto = 1e-6  # FIXME tendência de causar instabilidade numérica
    Rbaixa = 1e-6  # Ω
    Ralta = 1e6  # Ω
end


function testar_sem_falha_um_segmento()
    R_emissor_shunt = fill(-1.0, nc1)
    R_falha_shunt = fill(-1.0, num_pares_falhas)
    R_falha_serie = fill(-1.0, nc1)
    segmento_falha = 0
    comprimentos = [Lc]
    num_segmentos = length(comprimentos)
    aterrar_receptor = [true]
    sk, _ = laplace(ones(nt), tmax, nt)
    nf = length(sk)
    yn_segmentos = zeros(ComplexF64, (nc2, nc2, nf, num_segmentos))
    for k = 1:num_segmentos
        L1 = comprimentos[k]
        nome = "yn_" * string(round(L1, sigdigits=2))
        yn_segmentos[:,:,:,k] .= reshape(h5read(arquivo_matrizes, nome), (nc2, nc2, nf))
    end
    terminal_fonte = 1
    vout_t = simular_cabo(yn_segmentos, R_curto, R_fonte, R_emissor_shunt,
                          R_falha_shunt, R_falha_serie, segmento_falha,
                          pares_falhas, terminal_fonte, fases, blindagens,
                          armadura, v_fonte_t, tmax, nt, aterrar_receptor)
    arquivo_resultados = "teste/sem_falha_t$(terminal_fonte)_umseg.csv"
    salvar_resultados(arquivo_resultados, vout_t, v_fonte_t, tempo,
                      comprimentos, segmento_falha, pares_falhas, terminal_fonte,
                      R_emissor_shunt, R_falha_serie, R_falha_shunt, R_fonte,
                      nt_trunc, nome_condutores, aterrar_receptor)
    # A simetria geométrica faz com que a reposta dos condutores não excitados
    # sejam iguais.
    tol = 1e-3
    arquivo_resultados = "teste/sem_falha_t$(fases[1])_umseg.csv"
    df1 = CSV.read(arquivo_resultados, DataFrame; ntasks=1, silencewarnings=true)
    vf1 = df1[!,"v_" * nome_condutores[fases[1]]]
    vb1 = df1[!,"v_" * nome_condutores[blindagens[1]]]
    @assert any(sum(abs.(vf1)) > 0.0)
    vf2 = df1[!,"v_" * nome_condutores[fases[2]]]
    vf3 = df1[!,"v_" * nome_condutores[fases[3]]]
    @assert maximum(abs.(vf2 - vf3)) < tol
    vb2 = df1[!,"v_" * nome_condutores[blindagens[2]]]
    vb3 = df1[!,"v_" * nome_condutores[blindagens[3]]]
    @assert maximum(abs.(vb2 - vb3)) < tol
    return true
end
@assert testar_sem_falha_um_segmento()


function testar_sem_falha_dois_segmentos()
    R_emissor_shunt = fill(-1.0, nc1)
    R_falha_shunt = fill(-1.0, num_pares_falhas)
    R_falha_serie = fill(-1.0, nc1)
    segmento_falha = 0
    comprimentos = [Lc/2, Lc/2]
    num_segmentos = length(comprimentos)
    aterrar_receptor = [false, true]
    sk, _ = laplace(ones(nt), tmax, nt)
    nf = length(sk)
    yn_segmentos = zeros(ComplexF64, (nc2, nc2, nf, num_segmentos))
    for k = 1:num_segmentos
        L1 = comprimentos[k]
        nome = "yn_" * string(round(L1, sigdigits=2))
        yn_segmentos[:,:,:,k] .= reshape(h5read(arquivo_matrizes, nome), (nc2, nc2, nf))
    end
    for terminal_fonte in fases
        vout_t = simular_cabo(yn_segmentos, R_curto, R_fonte, R_emissor_shunt,
                            R_falha_shunt, R_falha_serie, segmento_falha,
                            pares_falhas, terminal_fonte, fases, blindagens,
                            armadura, v_fonte_t, tmax, nt, aterrar_receptor)
        arquivo_resultados = "teste/sem_falha_t$(terminal_fonte).csv"
        salvar_resultados(arquivo_resultados, vout_t, v_fonte_t, tempo,
                          comprimentos, segmento_falha, pares_falhas, terminal_fonte,
                          R_emissor_shunt, R_falha_serie, R_falha_shunt, R_fonte,
                          nt_trunc, nome_condutores, aterrar_receptor)
    end
    # tem que dar um resultado igual ao sem falha com um segmento
    tol = 1e-3
    arquivo_resultados = "teste/sem_falha_t$(fases[1])_umseg.csv"
    df1 = CSV.read(arquivo_resultados, DataFrame; ntasks=1, silencewarnings=true)
    vf1 = df1[!,"v_" * nome_condutores[fases[1]]]
    vb1 = df1[!,"v_" * nome_condutores[blindagens[1]]]
    for i in eachindex(fases)
        arquivo_resultados = "teste/sem_falha_t$(fases[i]).csv"
        df2 = CSV.read(arquivo_resultados, DataFrame; ntasks=1, silencewarnings=true)
        vf2 = df2[!,"v_" * nome_condutores[fases[i]]]
        vb2 = df2[!,"v_" * nome_condutores[blindagens[i]]]
        @assert maximum(abs.(vf1 - vf2)) < tol
        @assert maximum(abs.(vb1 - vb2)) < tol
    end
    return true
end
@assert testar_sem_falha_dois_segmentos()


function testar_falha_shunt_Ralta()
    R_emissor_shunt = fill(-1.0, nc1)
    R_falha_shunt = fill(-1.0, num_pares_falhas)
    R_falha_serie = fill(-1.0, nc1)
    segmento_falha = 1
    comprimentos = [Lc/2, Lc/2]
    num_segmentos = length(comprimentos)
    aterrar_receptor = [false, true]
    sk, _ = laplace(ones(nt), tmax, nt)
    nf = length(sk)
    yn_segmentos = zeros(ComplexF64, (nc2, nc2, nf, num_segmentos))
    for k = 1:num_segmentos
        L1 = comprimentos[k]
        nome = "yn_" * string(round(L1, sigdigits=2))
        yn_segmentos[:,:,:,k] .= reshape(h5read(arquivo_matrizes, nome), (nc2, nc2, nf))
    end
    for terminal_fonte in fases
        R_falha_shunt[:] .= -1.0
        R_falha_shunt[terminal_fonte] = Ralta
        vout_t = simular_cabo(yn_segmentos, R_curto, R_fonte, R_emissor_shunt,
                              R_falha_shunt, R_falha_serie, segmento_falha,
                              pares_falhas, terminal_fonte, fases, blindagens,
                              armadura, v_fonte_t, tmax, nt, aterrar_receptor)
        arquivo_resultados = "teste/falha_shunt_Ralta_t$(terminal_fonte).csv"
        salvar_resultados(arquivo_resultados, vout_t, v_fonte_t, tempo,
                          comprimentos, segmento_falha, pares_falhas, terminal_fonte,
                          R_emissor_shunt, R_falha_serie, R_falha_shunt, R_fonte,
                          nt_trunc, nome_condutores, aterrar_receptor)
    end
    # tem que dar um resultado similar ao sem falha
    tol = 1e-3
    arquivo_resultados = "teste/sem_falha_t$(fases[1])_umseg.csv"
    df1 = CSV.read(arquivo_resultados, DataFrame; ntasks=1, silencewarnings=true)
    vf1 = df1[!,"v_" * nome_condutores[fases[1]]]
    vb1 = df1[!,"v_" * nome_condutores[blindagens[1]]]
    for i in eachindex(fases)
        arquivo_resultados = "teste/falha_shunt_Ralta_t$(fases[i]).csv"
        df2 = CSV.read(arquivo_resultados, DataFrame; ntasks=1, silencewarnings=true)
        vf2 = df2[!,"v_" * nome_condutores[fases[i]]]
        vb2 = df2[!,"v_" * nome_condutores[blindagens[i]]]
        @assert maximum(abs.(vf1 - vf2)) < tol
        @assert maximum(abs.(vb1 - vb2)) < tol
    end
    return true
end
@assert testar_falha_shunt_Ralta()


function testar_falha_serie_Rbaixa()
    R_emissor_shunt = fill(-1.0, nc1)
    R_falha_shunt = fill(-1.0, num_pares_falhas)
    R_falha_serie = fill(-1.0, nc1)
    segmento_falha = 1
    comprimentos = [Lc/2, Lc/2]
    num_segmentos = length(comprimentos)
    aterrar_receptor = [false, true]
    sk, _ = laplace(ones(nt), tmax, nt)
    nf = length(sk)
    yn_segmentos = zeros(ComplexF64, (nc2, nc2, nf, num_segmentos))
    for k = 1:num_segmentos
        L1 = comprimentos[k]
        nome = "yn_" * string(round(L1, sigdigits=2))
        yn_segmentos[:,:,:,k] .= reshape(h5read(arquivo_matrizes, nome), (nc2, nc2, nf))
    end
    for terminal_fonte in fases
        R_falha_serie[:] .= -1.0
        R_falha_serie[terminal_fonte] = Rbaixa
        vout_t = simular_cabo(yn_segmentos, R_curto, R_fonte, R_emissor_shunt,
                              R_falha_shunt, R_falha_serie, segmento_falha,
                              pares_falhas, terminal_fonte, fases, blindagens,
                              armadura, v_fonte_t, tmax, nt, aterrar_receptor)
        arquivo_resultados = "teste/falha_serie_Rbaixa_t$(terminal_fonte).csv"
        salvar_resultados(arquivo_resultados, vout_t, v_fonte_t, tempo,
                          comprimentos, segmento_falha, pares_falhas, terminal_fonte,
                          R_emissor_shunt, R_falha_serie, R_falha_shunt, R_fonte,
                          nt_trunc, nome_condutores, aterrar_receptor)
    end
    # tem que dar um resultado similar ao sem falha
    tol = 1e-3
    arquivo_resultados = "teste/sem_falha_t$(fases[1])_umseg.csv"
    df1 = CSV.read(arquivo_resultados, DataFrame; ntasks=1, silencewarnings=true)
    vf1 = df1[!,"v_" * nome_condutores[fases[1]]]
    vb1 = df1[!,"v_" * nome_condutores[blindagens[1]]]
    for i in eachindex(fases)
        arquivo_resultados = "teste/falha_serie_Rbaixa_t$(fases[i]).csv"
        df2 = CSV.read(arquivo_resultados, DataFrame; ntasks=1, silencewarnings=true)
        vf2 = df2[!,"v_" * nome_condutores[fases[i]]]
        vb2 = df2[!,"v_" * nome_condutores[blindagens[i]]]
        @assert sqrt(sum((vf1 .- vf2).^2) / nt) < tol
        @assert sqrt(sum((vb1 .- vb2).^2) / nt) < tol
    end
    return true
end
@assert testar_falha_serie_Rbaixa()


function testar_ambas_falhas()
    R_emissor_shunt = fill(-1.0, nc1)
    R_falha_shunt = fill(-1.0, num_pares_falhas)
    R_falha_serie = fill(-1.0, nc1)
    segmento_falha = 1
    comprimentos = [Lc/2, Lc/2]
    num_segmentos = length(comprimentos)
    aterrar_receptor = [false, true]
    sk, _ = laplace(ones(nt), tmax, nt)
    nf = length(sk)
    yn_segmentos = zeros(ComplexF64, (nc2, nc2, nf, num_segmentos))
    for k = 1:num_segmentos
        L1 = comprimentos[k]
        nome = "yn_" * string(round(L1, sigdigits=2))
        yn_segmentos[:,:,:,k] .= reshape(h5read(arquivo_matrizes, nome), (nc2, nc2, nf))
    end
    for terminal_fonte in fases
        R_falha_shunt[:] .= -1.0
        R_falha_shunt[terminal_fonte] = Ralta
        R_falha_serie[:] .= -1.0
        R_falha_serie[terminal_fonte] = Rbaixa
        vout_t = simular_cabo(yn_segmentos, R_curto, R_fonte, R_emissor_shunt,
                              R_falha_shunt, R_falha_serie, segmento_falha,
                              pares_falhas, terminal_fonte, fases, blindagens,
                              armadura, v_fonte_t, tmax, nt, aterrar_receptor)
        arquivo_resultados = "teste/falha_ambas_falhas_t$(terminal_fonte).csv"
        salvar_resultados(arquivo_resultados, vout_t, v_fonte_t, tempo,
                          comprimentos, segmento_falha, pares_falhas, terminal_fonte,
                          R_emissor_shunt, R_falha_serie, R_falha_shunt, R_fonte,
                          nt_trunc, nome_condutores, aterrar_receptor)
    end
    # tem que dar um resultado similar ao sem falha
    tol = 1e-3
    arquivo_resultados = "teste/sem_falha_t$(fases[1])_umseg.csv"
    df1 = CSV.read(arquivo_resultados, DataFrame; ntasks=1, silencewarnings=true)
    vf1 = df1[!,"v_" * nome_condutores[fases[1]]]
    vb1 = df1[!,"v_" * nome_condutores[blindagens[1]]]
    for i in eachindex(fases)
        arquivo_resultados = "teste/falha_ambas_falhas_t$(fases[i]).csv"
        df2 = CSV.read(arquivo_resultados, DataFrame; ntasks=1, silencewarnings=true)
        vf2 = df2[!,"v_" * nome_condutores[fases[i]]]
        vb2 = df2[!,"v_" * nome_condutores[blindagens[i]]]
        @assert sqrt(sum((vf1 .- vf2).^2) / nt) < tol
        @assert sqrt(sum((vb1 .- vb2).^2) / nt) < tol
    end
    return true
end
@assert testar_ambas_falhas()

# deletar resultados
foreach(rm, filter(endswith(".csv"), readdir("teste",join=true)))
rm("teste/umbilical_tripolar.h5")
println("Testes da simulação do cabo bem sucedidos.")


#=
using Plots

function lerplotar_resultados(arquivo_resultados)
    df = CSV.read(arquivo_resultados, DataFrame)
    t = df[!,"tempo"] * 1e3
    v_fase_1 = df[!,"v_fase_1"]
    v_blindagem_1 = df[!,"v_blindagem_1"]
    v_fase_2 = df[!,"v_fase_2"]
    v_blindagem_2 = df[!,"v_blindagem_2"]
    v_fase_3 = df[!,"v_fase_3"]
    v_blindagem_3 = df[!,"v_blindagem_3"]
    v_armadura = df[!,"v_armadura"]
    p = plot(xlabel="Tempo [ms]", ylabel="Tensão [V]")
    plot!(t, v_fase_1, label="fase 1", color=:blue, linestyle=:solid)
    plot!(t, v_blindagem_1, label="blindagem 1", color=:blue, linestyle=:dash)
    plot!(t, v_fase_2, label="fase 2", color=:darkorange, linestyle=:solid)
    plot!(t, v_blindagem_2, label="blindagem 2", color=:darkorange, linestyle=:dash)
    plot!(t, v_fase_3, label="fase 3", color=:darkgreen, linestyle=:solid)
    plot!(t, v_blindagem_3, label="blindagem 3", color=:darkgreen, linestyle=:dash)
    plot!(t, v_armadura, label="armadura", color=:black, linestyle=:dashdot)
    return p
end

lerplotar_resultados(arquivo_resultados)
=#