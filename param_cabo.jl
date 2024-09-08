#=
Calcula os parâmetros elétricos de um cabo umbilical.
=#
using SpecialFunctions
using LinearAlgebra
using HDF5
using CSV, DataFrames

const MU_0 = 4e-7 * pi  # permeabilidade magnética do vácuo
const EPSILON_0 = 8.8541878128e-12  # permissividade elétrica

# %% Parâmetros dos cabos

"""Matrizes Z e Y do cabo com blindagem e isolação.

Parâmetros
----------
    w : frequência angular [rad/s].
    rc : vetor de raios do cabo [m]: condutor central, isolante, blindagem e isolante externo.
    rho_c : resistividade elétrica do condutor central [Ohm.m].
    rho_b : resistividade elétrica da blindagem [Ohm.m].
    epsr_c : permissividade elétrica relativa do condutor central -- corrigida para
             levar em conta o efeito da camada semicondutora.
    epsr_b : permissividade elétrica relativa da capa externa, jaqueta, da blindagem da veia
    mur_c : permeabilidade magnética relativa do condutor central.
    mur_b : permeabilidade magnética relativa da capa externa, jaqueta, da blindagem da veia.
    td1 : tangente de perdas do dielétrico 1
    td2 : tangente de perdas do dielétrico 2

Retorna
-------
    Z : matriz de impedância longitudinal [Ohm/m]
    Y : matriz de admitância transversal [S/m]
"""
function cZYcbi(w, rc, rho_c, rho_b, epsr_c, epsr_b, mur_c=1, mur_b=1, td1=0.0001, td2=0.0004)
    jwu0 = 1.0im * w * MU_0
    jwu0_2pi = jwu0 / (2 * pi)
    etac = sqrt((jwu0 * mur_c) / rho_c)
    etab = sqrt((jwu0 * mur_b) / rho_b)
    z1 = etac * rho_c / (2 * pi * rc[1]) * besselix(0, etac * rc[1]) / besselix(1, etac * rc[1])
    z2 = jwu0_2pi * log(rc[2] / rc[1])
    eta1 = etab * rc[2]
    eta2 = etab * rc[3]
    escala = exp(abs(real(eta1)) - eta2 - abs(real(eta2)) + eta1)  # s1/s2
    den = (besselix(1, eta2) * besselkx(1, eta1)) - (besselix(1, eta1) * besselkx(1, eta2)) * escala
    num1 = (besselix(0, eta1) * besselkx(1, eta2)) * escala + (besselix(1, eta2) * besselkx(0, eta1))
    z3 = (rho_b * etab) / (2 * pi * rc[2]) * num1 / den
    s2 = exp(abs(real(eta2)) - eta1)
    if isinf(s2)
        z4 = 0.0
    else
        z4 = rho_b / (2 * pi * rc[2] * rc[3] * den * s2)
    end
    num2 = (besselix(0, eta2) * besselkx(1, eta1)) + (besselix(1, eta1) * besselkx(0, eta2)) * escala
    z5 = (rho_b * etab) / (2 * pi * rc[3]) * num2 / den
    z6 = jwu0_2pi * log(rc[4] / rc[3])
    Z = zeros(ComplexF64, (2,2))
    Z[1,1] = z1 + z2 + z3 + z5 + z6 - 2 * z4
    Z[1,2] = z5 + z6 - z4
    Z[2,1] = z5 + z6 - z4
    Z[2,2] = z5 + z6
    # g não depende da frequência pois o valor de td é obtido por ensaios em uma frequência
    y1 = 1.0im * w * 2 * pi * (epsr_c * EPSILON_0) / log(rc[2] / rc[1])
    y2 = 1.0im * w * 2 * pi * (epsr_b * EPSILON_0) / log(rc[4] / rc[3])
    g1 = td1 * 120 * pi^2 * (epsr_c * EPSILON_0) / log(rc[2] / rc[1])
    g2 = td1 * 120 * pi^2 * (epsr_c * EPSILON_0) / log(rc[2] / rc[1])
    y1 = 1.0im * w * 2 * pi * (epsr_c * EPSILON_0) / log(rc[2] / rc[1]) + g1
    y2 = 1.0im * w * 2 * pi * (epsr_b * EPSILON_0) / log(rc[4] / rc[3]) + g2
    Y = zeros(ComplexF64, (2,2))
    Y[1,1] = y1
    Y[1,2] = -y1
    Y[2,1] = -y1
    Y[2,2] = y1 + y2
    return Z, Y
end


"""Matrizes Z e P próprias do cabo Pipe-Type (PT) em função da camada metálica.

Parâmetros
----------
    w : frequência angular [rad/s].
    mur_a : permeabilidade magnética relativa da armadura.
    rho_a : resistividade elétrica da armadura [Ohm.m].
    epsr_ae : permissividade elétrica da última camada isolante, que envolve a armadura.
    theta_jk : ângulo entre os cabos coaxiais [rad].
    rc : vetor de raios do cabo [m]: condutor central, isolante, blindagem e isolante externo.
    ra : vetor de raios da armadura [m]: interno, externo, capa.
    di : vetor de distâncias do centro do umbilical ao centro das veias [m].

Retorna
-------
    Z : matriz de impedância longitudinal [Ohm/m]
    P : matriz de coeficientes potenciais de Maxwell [V/c/m]
"""
function cZPpipein(w, mur_a, rho_a, epsr_ae, theta_jk, rc, ra, di)
    ncx = 2
    nn = 24
    re = rc[4]
    rp1 = ra[1]
    y1 = rp1 * sqrt((1im * MU_0 * mur_a * w) / rho_a)
    q = zeros(ComplexF64, (length(di), length(di)))
    for linha = 1:length(di)
        for coluna = 1:length(di)
            sum = 0.0
            sum2 = 0.0
            if linha == coluna
                for n = 1:nn
                    sum = sum + (di[linha] / rp1)^(2*n) / ((mur_a + 1) * n + y1 * besselkx(n-1, y1) / besselkx(n, y1))
                end
                sum = 2.0 * mur_a * sum + log(0.0im + (rp1 * (1 - (di[linha] / rp1)^2)) / re)
                q[linha,coluna] = sum
                sum = 0.0
            else
                for n = 1:nn
                    sum = sum - ((cos(theta_jk * n) * ((di[linha] * di[coluna]) / rp1^2)^n) / n)
                    sum2 = sum2 + (cos(theta_jk * n) * ((di[linha] * di[coluna]) / rp1^2)^n) / ((mur_a + 1) * n + y1 * besselkx(n-1, y1) / besselkx(n, y1))
                end
                q[linha,coluna] = sum + 2*mur_a*sum2 + log(0.0im + rp1 / sqrt(-(2 * cos(theta_jk) * di[linha] * di[coluna]) + di[linha]^2 + di[coluna]^2))
                sum = 0.0
                sum2 = 0.0
            end
        end
    end
    zpp = @. ((1im * MU_0 * w) * (q + mur_a / y1 * besselkx(0, y1) / besselkx(1, y1))) / (2 * pi)
    auxZ = zeros(ComplexF64, size(zpp) .* ncx)
    for linha = 0:(size(zpp)[1] - 1)
        i1 = (linha) * ncx + 1
        i2 = i1 + ncx - 1
        for coluna = 0:(size(zpp)[2] - 1)
            k1 = (coluna) * ncx + 1
            k2 = k1 + ncx - 1
            auxZ[i1:i2, k1:k2] .= zpp[linha+1, coluna+1]
        end
    end
    Zpipe = zeros(ComplexF64, size(auxZ) .+ 1)
    for linha = 1:size(auxZ)[1]
        for coluna = 1:size(auxZ)[2]
            Zpipe[linha,coluna] = auxZ[linha,coluna]
        end
    end
    qp = zeros(ComplexF64, (length(di),length(di)))
    for linha = 1:length(di)
        for coluna = 1:length(di)
            sum = 0.0
            sum2 = 0.0
            if linha==coluna
                qp[linha,coluna] = log(0.0im + (rp1 * (1 - (di[linha] / rp1)^2)) / re)
            else
                for n = 1:nn
                    sum = sum - ((cos(theta_jk * n) * ((di[linha] * di[coluna]) / rp1^2)^n) / n)
                end
                qp[linha,coluna] = sum + log(0.0im + rp1 / sqrt(-(2 * cos(theta_jk) * di[linha] * di[coluna]) + di[linha]^2 + di[coluna]^2))
                sum = 0.0
                sum2 = 0.0
            end
        end          
    end
    ppp = qp / (2 * pi * EPSILON_0 * epsr_ae)
    auxP = zeros(ComplexF64, size(ppp) .* ncx)
    for linha = 0:(size(ppp)[1] - 1)
        i1 = (linha) * ncx + 1
        i2 = i1 + ncx - 1
        for coluna = 0:(size(ppp)[2] - 1)
            k1 = (coluna) * ncx + 1
            k2 = k1 + ncx - 1
            auxP[i1:i2, k1:k2] .= ppp[linha+1, coluna+1]
        end
    end
    Ppipe = zeros(ComplexF64, size(auxP) .+ 1)
    for linha = 1:size(auxP)[1]
        for coluna = 1:size(auxP)[2]
            Ppipe[linha,coluna] = auxP[linha,coluna]
        end
    end
    return Zpipe, Ppipe
end


"""Matrizes Z e P cabo própria da armadura Pipe-Type (PT).

Parâmetros
----------
    w : frequência angular [rad/s].
    ra : vetor de raios da armadura [m]: interno, externo, capa.
    mur_a : permeabilidade magnética relativa da armadura.
    rho_a : resistividade elétrica da armadura [Ohm.m].
    epsr_ae : permissividade elétrica da última camada isolante, que envolve a armadura.

Retorna
-------
    Z : matriz de impedância longitudinal [Ohm/m]
    P : matriz de coeficientes potenciais de Maxwell [V/c/m]
"""
function cZPpipe(w, ra, mur_a, rho_a, epsr_ae)
    ncx = 2  # número de condutores "ativos" em cada cabo SC (single core)
    rae = ra[2]
    rai = ra[1]
    ma = sqrt((1.0im * w * MU_0 * mur_a) / rho_a)
    eta_i = ma * rai
    eta_e = ma * rae
    escala = exp(abs(real(eta_i)) - eta_e - abs(real(eta_e)) + eta_i)  # si/se
    den = (besselix(1, eta_e) * besselkx(1, eta_i)) - (besselix(1, eta_i) * besselkx(1, eta_e)) * escala
    numi = (besselix(0, eta_i) * besselkx(1, eta_e)) * escala + (besselix(1, eta_e) * besselkx(0, eta_i))
    nume = (besselix(0, eta_e) * besselkx(1, eta_i)) + (besselix(1, eta_i) * besselkx(0, eta_e)) * escala
    zai = (ma * rho_a) / (2 * pi * rai) * numi / den
    zae = (ma * rho_a) / (2 * pi * rae) * nume / den
    se = exp(abs(real(eta_e)) - eta_i)
    if isinf(se)
        zam = 0.0
    else
        zam = rho_a / (2 * pi * rai * rae * den * se)
    end
    rf = ra[3]
    zins = ((1.0im * w * MU_0) / (2 * pi)) * log(rf / rae)
    zc1 = zai + zae - 2*zam + zins
    zc2 = zae - zam + zins
    zc3 = zae + zins
    Zp = zc1
    Zpipe = zeros(ComplexF64, (7,7))  # FIXME e quando houver mais de 3 fases?
    for i = 1:length(ra)
        for j = 1:length(ra)
            Zpipe[2*i-1:2*i, 2*j-1:2*j] .= Zp
        end
    end
    Zpipe[7, 1:6] .= zc2
    Zpipe[1:6, 7] .= zc2
    Zpipe[7, 7] = zc3
    # Zpipe= [Zp Zp Zp zc2 Zp Zp Zp zc2 Zp Zp Zp zc2 zc2 zc2 zc2 zc3]
    pins = log(rf / rae) / (2 * pi * EPSILON_0 * epsr_ae)
    n = 3 * ncx + 1
    Ppipe = pins * ones((n, n))
    return Zpipe, Ppipe
end


"""Impedância de retorno pelo mar do cabo pipe-Type (PT).

Parâmetros
----------
    w : frequência angular [rad/s].
    rca : raio da capa externa da armadura [m].
    sigma : permissividade elétrica do mar [S/m].
    epsr : permissividade elétrica do mar.

Retorna
-------
    Z : matriz de impedância longitudinal [Ohm/m]
"""
function cZmar(w, rca, sigma=5.0, epsr=81.0)
    rho = 1.0 / sigma
    eta = sqrt(1.0im * w * MU_0 * (sigma + 1.0im * w * EPSILON_0 * epsr))
    Zmar = eta * rho / (2 * pi * rca) * besselkx(0, eta * rca) / besselkx(1, eta * rca)
    return Zmar
end


"""Elimina os cabos pára-raios (ou blindagem) da matriz, considerando-os
continuamente aterrados (redução de Kron).

Parâmetros
----------
    M : matriz na qual fazer a redução de Kron.
    np : número de cabos a eliminar, considerados os últimos da matriz.

Retorna
-------
    Msprc : matriz reduzida
"""
function eliprc(M, np)
    Minv = inv(M)
    lenM = size(M)[1]
    Msprc = Minv[1:lenM-np, 1:lenM-np]
    return Msprc
end


"""Calcula a matriz de admitância nodal Yn.
    Yn = [[y11; y12]
            [y12; y11]]

Parâmetros
----------
    zc : matriz de impedância longitudinal [Ohm/m].
    yc : matriz de admitância transversal [S/m].
    L : comprimento total [m].

Retorna
-------
    yn : [[y11; y12], [y12, y11]]
"""
function ynodal(zc, yc, L)
    Z = zc
    Y = yc
    n = size(zc)[1]
    sqrtzy = sqrt(Z * Y)
    Yc = inv(Z) * sqrtzy
    H = exp(-L * sqrtzy)
    H2 = H * H
    I = diagm(ones(n))
    inv_IH2 = inv(I - H2)
    YL1 = Yc * (I + H2) * inv_IH2
    YL2 = -2.0 * Yc * H * inv_IH2
    yn = [[YL1  YL2]; [YL2  YL1]]
    yn = (yn + transpose(yn)) / 2.0
    #@assert issymmetric(yn)
    return yn
end


"""Cacula as matrizes de impedância Z e admitância Y por unidade de
comprimento de um cabo umbilical.

Parâmetros
----------
    w : frequência angular [rad/s].
    rc : vetor de raios do cabo [m]: condutor central, isolante, blindagem e isolante externo.
    ra : vetor de raios da armadura [m]: interno, externo, capa.
    rho_a : resistividade elétrica da armadura [Ohm.m].
    rho_b : resistividade elétrica da blindagem [Ohm.m].
    rho_c : resistividade elétrica do condutor central [Ohm.m].
    epsr_a : permissividade elétrica da última camada isolante, que envolve a armadura.
    epsr_b : permissividade elétrica relativa da capa externa, jaqueta, da blindagem da veia
    epsr_c : permissividade elétrica relativa do condutor central -- corrigida para
             levar em conta o efeito da camada semicondutora.
    mur_a : permeabilidade magnética relativa da armadura.
    theta_jk : ângulo entre os cabos coaxiais [rad].
    di : vetor de distâncias do centro do umbilical ao centro das veias [m].
    sig_s : condutividade elétrica do mar [S/m].
    eps_s : permissividade elétrica relativa do mar.
    nx : número de condutores "ativos" em cada cabo SC (single core)

Retorna
-------
    Zcabo : matriz de impedância série por unidade de comprimento [Ohm/m].
    Ycabo : matriz de admitância shunt por unidade de comprimento [S/m].

"""
function zy_cabo(w, rc, ra, rho_a, rho_b, rho_c, epsr_a, epsr_b, epsr_c, mur_a,
                 theta_jk, di, sig_s, eps_s, nx=2)
    Zin, Yin = cZYcbi(w, rc, rho_c, rho_b, epsr_c, epsr_b, 1, 1)
    Zpipein, Ppipein = cZPpipein(w, mur_a, rho_a, epsr_a, theta_jk, rc, ra, di)
    Zpipe, Ppipe = cZPpipe(w, ra, mur_a, rho_a, epsr_a)
    # Calculo de Z
    Z0 = cZmar(w, ra[end], sig_s, eps_s) * ones((3*nx + 1, 3*nx + 1))
    Zint = zeros(ComplexF64, size(Z0))
    for i = 1:3
        Zint[2*i-1:2*i, 2*i-1:2*i] .= Zin
    end
    Zcabo = Zint + Zpipein + Zpipe + Z0
    
    pin = inv(Yin) * (1.0im * w)
    pint = zeros(ComplexF64, (size(Ppipe)[1], size(Ppipe)[1]))
    for i = 1:3
        pint[2*i-1:2*i, 2*i-1:2*i] .= pin
    end
    Ycabo = 1.0im * w * inv(pint + Ppipein + Ppipe)
    # Ycabo não está simétrica, mas `maximum(abs.(Ycabo - transpose(Ycabo)))`
    # é da ordem de 1e-20. Isto é resultado da inversão, pois
    # `issymmetric(pint + Ppipein + Ppipe) == true`.
    Ycabo = (Ycabo + transpose(Ycabo)) / 2.0  # forçar simetria
    #@assert issymmetric(Zcabo)
    #@assert issymmetric(Ycabo)
    return Zcabo, Ycabo
end


"""Rotacao da Matriz de Autovetores para minimizar parte Imaginaria. 
Complemento ao Metodo de Newton-Raphson

Parâmetros
----------
    S : Matriz de autovetores.

Retorna
-------
    rot(S) : matriz com parte imaginária minimizada.

"""
function rot(S)
    Nc = size(S)[1]
    scale1 = zeros(ComplexF64, Nc)
    scale2 = copy(scale1)
    scale = copy(scale1)
    ang = copy(scale1)
    err1 = copy(scale1)
    err2 = copy(scale1)
    numerator = zeros(Nc)
    denominator = copy(numerator)
    SA = zeros(ComplexF64, (Nc, Nc))
    SB = copy(SA)
    for col = 1:Nc
        for j = 1:Nc
            numerator[col] += imag(S[j, col]) * real(S[j, col])
            denominator[col] += real(S[j, col])^2 - imag(S[j, col])^2
        end
        numerator[col] *= -2
        ang[col] = 0.5 * atan(denominator[col], numerator[col])
        scale1[col] = cos(ang[col]) + 1.0im * sin(ang[col])
        scale2[col] = cos(ang[col] + pi / 2) + 1.0im * sin(ang[col] + pi / 2)
        aaa = bbb = ccc = ddd = eee = fff = 0.0
        
        for j = 1:Nc
            SA[j, col] *= scale1[col]
            SB[j, col] *= scale2[col]
            aaa += imag(SA[j, col])^2
            bbb += real(SA[j, col]) * imag(SA[j, col])
            ccc += real(SA[j, col])^2
            ddd += imag(SB[j, col])^2
            eee += real(SB[j, col]) * imag(SB[j, col])
            fff += real(SB[j, col])^2
        end
        err1[col] = aaa * cos(ang[col])^2 + bbb * sin(2 * ang[col]) + ccc * sin(ang[col])^2
        err2[col] = ddd * cos(ang[col])^2 + eee * sin(2 * ang[col]) + fff * sin(ang[col])^2
        if abs(err1[col]) < abs(err2[col])
            scale[col] = scale1[col]
        else
            scale[col] = scale2[col]
        end
        S[:, col] *= scale[col]
    end
    return S
end


"""Elimina o Cruzamento de Autovetores. Método de Newton-Raphson.

Parâmetros
----------
    A : matriz
    V0 : autovetores
    d0 : autovalores
    tol : toleância para considerar convergência
    maxiter : número máximo de iterações

Retorna
-------
    outd :
    outv : 
"""
function eigNR(A, V0, d0, tol=1e-9, maxiter=10000)
    nc = size(A)[1]
    scale = norm(A)
    A /= scale
    outd = zeros(ComplexF64, nc)
    outv = zeros(ComplexF64, (nc, nc))
    jac = zeros(ComplexF64, (nc + 1, nc + 1))
    vec1 = zeros(ComplexF64, nc+1)
    res1 = zeros(ComplexF64, nc+1)
    for i = 1:nc
        v = V0[:,i]
        d = d0[i]
        resi = maximum(abs.(A * v - v .* d))
        itr = 0
        while resi > tol && itr < maxiter
            itr += 1
            jac[1:nc, 1:nc] = A - d * I
            jac[1:nc, (nc+1)] = -v
            jac[(nc+1), 1:nc] = 2 * v
            # dot(v,w) em julia conjuga v
            res1[:] = [A * v - v .* d ; dot(conj(v), v) - 1.0]
            for i = 1:(nc+1)
                jac[i,i] += 1e-12  # pra evitar divisão por zero
            end
            vec1[:] = qr(jac) \ res1
            resi = maximum(abs.(vec1))
            v = v - vec1[1:nc]
            d = d - vec1[nc+1]
        end
        if itr == maxiter
            println("máximo de iterações em eigNR")
        end
        outd[i] = d * scale
        outv[:, i] = v ./ norm(v)
    end
    return outd, outv
end
