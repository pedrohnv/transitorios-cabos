#=
Transformada numérica de Laplace.

Referências
-----------
    GÓMEZ, Pablo; URIBE, Felipe A. The numerical Laplace transform: An accurate
    technique for analyzing electromagnetic transients on power system devices.
    International Journal of Electrical Power & Energy Systems, v. 31, n. 2-3,
    p. 116-123, 2009.
=#

using FFTW

"""Transformada de Laplace do vetor `y(t)`.

Parâmetros
----------
    y : sinal a ser transformado, assumido ser números Reais.
    tmax : último tempo [s].
    nt : número de passos no tempo.

Retorna
-------
    L(y) : vetor transformado.
    s : vetor de frequências complexas `s = c + jω`.
"""
function laplace(y::Array{Float64}, tmax::Float64, nt::Int)
    c = log(nt^2) / tmax
    dt = tmax / (nt - 1)
    dw = 2pi / tmax
    ns = (nt ÷ 2) + 1
    s = [c + 1im * dw * (k - 1) for k = 1:ns]
    v = [dt * exp(-c * (k - 1) * dt) * y[k] for k = 1:nt]
    return rfft(v), s
end


"""Transformada inversa de Laplace do vetor `y(s)`.

Parâmetros
----------
    y : sinal a ser transformado, assumido como pares de números Complexos conjugados.
    s : vetor de frequências complexas [rad/s].
    tmax : último tempo [s].
    nt : número de passos no tempo.
    filtro : Veja Tabela 1 de GÓMEZ & URIBE (2009):
        "Blackman" : σ(ω) = 0.42 + 0.5 cos(π ω/Ω) + 0.08 cos(2π ω/Ω)
        "Hanning" : σ(ω) = (1 +  cos(π ω/Ω)) / 2
        "Lanczos" : σ(ω) = sin(π ω/Ω) / (π ω/Ω)
        "Riesz" : σ(ω) = 1 - |ω/Ω|^2
        outros : σ(ω) = 1

Retorna
-------
    L^(-1)(y) : vetor transformado.
"""
function invlaplace(y, s, tmax, nt::Int, filtro="")
    filtro = uppercase(filtro)
    omega = imag(s[end])
    function sigma(w)
        alpha = pi * w / omega
        if filtro == uppercase("Blackman")
            return 0.42 + 0.5 * cos(alpha) + 0.08 * cos(2 * alpha)
        elseif filtro == uppercase("Hanning")
            return (1 + cos(pi * w / omega)) / 2
        elseif filtro == uppercase("Lanczos")
            return sin(alpha) / alpha
        elseif filtro == uppercase("Riesz")
            return 1.0 - abs(w / omega)^2.0
        elseif filtro == ""
            return 1.0
        else
            throw(ArgumentError("filtro não reconhecido: $(filtro)"))
        end
    end
    ns = length(s)
    g = zeros(Complex{Float64}, ns)
    for k = 1:ns
        w = imag(s[k])
        g[k] = sigma(w) * y[k]
    end
    v = irfft(g, nt)
    c = real(s[1])
    dt = tmax / (nt - 1)
    f = zeros(nt)
    for k = 1:nt
        f[k] = v[k] * exp(c * (k - 1) * dt) / (dt)
    end
    return f
end