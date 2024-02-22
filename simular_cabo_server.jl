#= Simular um cabo umbilical. Parâmetros lidos de um arquivo CSV. =#
using HDF5
using Parameters
using CSV, DataFrames

include("param_cabo.jl");


begin
    ## Geometria do umbilical PT
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

    fases = [1, 3, 5]  # não alterar
    blindagens = [2, 4, 6]  # não alterar
    armadura = 7  # não alterar
    condutores = [fases; blindagens; armadura]
    nc1 = length(condutores)  # número de condutores
    nc2 = nc1 * 2
    nome_condutores = sort(Dict([
        [fases[i] => "fase_$(i)" for i in eachindex(fases)];
        [blindagens[i] => "blindagem_$(i)" for i in eachindex(blindagens)];
        [armadura => "armadura"]
    ]))
end


function ignorar_missing(coluna)
    n = 0
    for i in eachindex(param[!, coluna])
        try
            Float64(param[i, coluna])
            n += 1
        catch
            break
        end
    end
    return Float64.(param[1:n,coluna])
end


begin  # ler parâmetros
    param = CSV.read("parametros.csv", DataFrame)
    dt = Float64.(param[1,"dt"])
    nt_trunc = Int.(param[1,"nt"])
    nt = Int(ceil(nt_trunc / 0.80))  # 20% dos tempos finais é lixo numérico
    tempo = range(0, step=dt, length=nt)
    tmax = tempo[end]
    comprimentos = ignorar_missing("comprimentos")
    num_comprimentos = length(comprimentos)
    aterrar_receptor = Bool.(ignorar_missing("aterrar_receptor"))
    @assert num_comprimentos == length(aterrar_receptor)
    segmento_falha = Int.(param[1,"segmento_falha"])
    @assert segmento_falha <= num_comprimentos
    terminal_fonte = Int.(param[1,"terminal_fonte"])
    @assert 1 <= terminal_fonte <= 7
    tmax_trunc = tempo[nt_trunc]  # tempo [s] final que será salvo
    R_emissor_shunt = ignorar_missing("R_emissor_shunt")
    R_falha_serie = ignorar_missing("R_falha_serie")
    R_chave = Float64.(param[1,"R_chave"])
    R_curto = Float64.(param[1,"R_curto"])
    R_falha_shunt = ignorar_missing("R_falha_shunt")
    num_pares_falhas = length(R_falha_shunt)
    pares_falhas = Vector{Vector{Int64}}(undef, num_pares_falhas)
    for i = 1:num_pares_falhas
        pares_falhas[i] = [param[i, "par_falha_1"], param[i, "par_falha_2"]]
    end
    arquivo_resultados = param[1, "arquivo_resultados"]
end


begin  # Fonte de tensão [V] lida de um arquivo. Cada valor em uma linha.
    arquivo_fonte = "v_fonte.txt"
    v_fonte_file = split(read(arquivo_fonte, String), "\n")
    v_fonte_t = zeros(nt)
    local zero_pad = false
    for i = 1:nt
        try
            v_fonte_t[i] = parse(Float64, v_fonte_file[i])
        catch e 
            if e isa BoundsError || e isa ArgumentError
                v_fonte_t[i] = 0.0
                zero_pad = true
            end
        end
    end
end


# %% Cálculos
arquivo_matrizes = "umbilical_tripolar.h5"
lk = ReentrantLock()
@lock lk precalc_matrizes(arquivo_matrizes, comprimentos, tmax, nt, rc, ra, rho_a, rho_b, rho_c, epsr_a, epsr_b, epsr_c, mur_a, theta_jk, di, sig_s, eps_s, nx)

# Carregar matrizes YN
sk, v_fonte_s = laplace(v_fonte_t, tmax, nt)
nf = length(sk)  # número de frequências
yn_comprimentos = zeros(ComplexF64, (nc2, nc2, nf, num_comprimentos))
for k = 1:num_comprimentos
    L1 = comprimentos[k]
    nome = "yn_" * string(round(L1, sigdigits=2))
    @lock lk yn_comprimentos[:,:,:,k] .= reshape(h5read(arquivo_matrizes, nome), (nc2, nc2, nf))
end

vout_t = simular_cabo(yn_comprimentos, R_curto, R_chave, R_emissor_shunt, R_falha_shunt, R_falha_serie, segmento_falha, pares_falhas, terminal_fonte, fases, blindagens, armadura, v_fonte_t, tmax, nt, aterrar_receptor)
@lock lk salvar_resultados(arquivo_resultados, vout_t, v_fonte_t, tempo, comprimentos, segmento_falha, pares_falhas, terminal_fonte, R_emissor_shunt, R_falha_serie, R_falha_shunt, R_chave, nt_trunc, nome_condutores, aterrar_receptor)
