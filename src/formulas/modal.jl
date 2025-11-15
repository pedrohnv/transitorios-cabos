#= Functions for modal analysis. =#

using LinearAlgebra
using PythonCall

least_squares = pyimport("scipy.optimize").least_squares


"""Solve the complex symmetric eigenvalue problem for multiple frequencies.

S_omega may need to be normalized. Example:
```
w2u0e0 = -omega.^2 * mu_0 * epsilon_0
S_omega = [A[i] / (-w2u0e0[i]) - I for i in 1:K]

eigvals, eigvecs = eig_levenberg_marquardt(S_omega, tol=1e-11, max_iter=10000)

# scale back the eigenvalues; eigenvectors remain unchanged
eigvals = [(-w2u0e0 .* (1 .+ eigvals[:, i])) for i in 1:N]
```

# Parameters
- `S_omega`: (K x N x N) array of complex matrices S(ω) for K frequencies
- `tol`: Tolerance for convergence
- `max_iter`: Maximum number of iterations

# Returns
- `eigenvalues`: (K x N) array of complex eigenvalues
- `eigenvectors`: (K x N x N) array of complex eigenvectors (columns are eigenvectors)

# References

A. I. Chrysochos, T. A. Papadopoulos and G. K. Papagiannis, "Robust Calculation
of Frequency-Dependent Transmission-Line Transformation Matrices Using the
Levenberg–Marquardt Method," in IEEE Transactions on Power Delivery, vol. 29,
no. 4, pp. 1621-1629, Aug. 2014, doi: 10.1109/TPWRD.2013.2284504.
"""
function eig_levenberg_marquardt(S_omega; tol=1e-8, max_iter=1000)
    N, _, K = size(S_omega)
    eigenvalues = zeros(ComplexF64, N, K)
    eigenvectors = zeros(ComplexF64, N, N, K)

    # First solve the first frequency with standard eigensolver
    S0 = S_omega[:, :, 1]
    F = eigen(S0)
    eigenvalues[:, 1] = F.values
    eigenvectors[:, :, 1] = F.vectors

    # Prepare real-valued formulation for subsequent frequencies
    for k in 2:K
        S = S_omega[:, :, k]
        S_re = real(S)
        S_im = imag(S)

        # Previous solution as initial guess
        prev_eigvecs = eigenvectors[:, :, k-1]
        prev_eigvals = eigenvalues[:, k-1]

        # Solve for each eigenvalue/eigenvector pair
        for i in 1:N
            # Initial guess from previous frequency
            t_prev = prev_eigvecs[:, i]
            lambda_prev = prev_eigvals[i]

            # Real-valued residual function
            function residuals(x)
                t_re = x[1:N]
                t_im = x[N+1:2*N]
                lambda_re = x[2*N+1]
                lambda_im = x[2*N+2]

                # Eigen equation residuals
                res1 = (S_re * t_re - S_im * t_im) - (lambda_re .* t_re - lambda_im .* t_im)
                res2 = (S_im * t_re + S_re * t_im) - (lambda_im .* t_re + lambda_re .* t_im)

                # Normalization constraints
                norm1 = sum(t_re.^2) - sum(t_im.^2) - 1  # |t|^2 = 1
                norm2 = 2 * sum(t_re .* t_im)  # Ensures proper phase

                return vcat(res1, res2, [norm1, norm2])
            end

            # Initial guess
            x0 = vcat(
                real(t_prev),
                imag(t_prev),
                [real(lambda_prev), imag(lambda_prev)]
            )

            # Normalize initial guess to satisfy constraints
            t_re = x0[1:N]
            t_im = x0[N+1:2*N]
            denom = sqrt(sum(t_re.^2) + sum(t_im.^2))
            x0[1:2*N] ./= denom
            
            # Solve with Levenberg-Marquardt
            res = least_squares(
                @py residuals,
                x0,
                method="lm",
                xtol=tol,
                ftol=tol,
                max_nfev=max_iter,
            )

            # Extract solution
            if !Bool(res.success)
                msg = "Warning: Did not converge for frequency $(k), eigenvalue $(i)\nResidual norm: $(norm(res.fun))"
                warnings.warn(msg)
            end

            # Extract solution
            res_x = pyconvert(Array{ComplexF64}, res.x)
            t_re = res_x[1:N]
            t_im = res_x[N+1 : 2N]
            lambda_re = res_x[2N + 1]
            lambda_im = res_x[2N + 2]
            
            # Store solution
            eigenvectors[:, i, k] = t_re .+ 1im .* t_im
            eigenvalues[i, k] = lambda_re + 1im * lambda_im

            # Ensure eigenvectors are properly normalized
            eigvec = eigenvectors[:, i, k]
            eigenvectors[:, i, k] = eigvec / sqrt(sum(abs.(eigvec).^2))
        end
    end

    return eigenvalues, eigenvectors
end


"""Calculates the propagation modes, their velocity and attenuation.

The attenutation `\\alpha` \\[Np/m\\] is the `\\ln` proportion of a wave at the
start `v_0` by the wave 1 m forward `v_1`

```math
\\alpha = \\ln(v_0 / v_1)
```
```math
v_1 = v_0 / \\exp(\\alpha)
```

# Parameters

- `Z`: Series impedance matrix per unit length \\[Ω/m\\]. Dimensões (Nc, Nc, Nf).
- `Y`: Shunt admittance matrix per unit length \\[S/m\\]. Dimensões (Nc, Nc, Nf).
- `complex_frequencies`: Vetor de frequências angulares complexas no formato `c + jw` [rad/s]. Dimensão (Nf,).
- `unwrap`: Unwrap the propagation modes so their order is preserved?
- `tol`: Tolerance for convergence if unwrap is true.
- `max_iter`: Maximum number of iterations if unwrap is true.

# Returns

- `propagation`: Propagation modes `a + jb`, where `a` is \\[Np/m\\] and `b` is \\[rad/s\\]. Dimensões (Nc, Nf).
- `velocity`: Velocity of the propagation modes in \\[m/μs\\]. Dimensões (Nc, Nf).
- `attenuation`: Attenuation of the propagation modes in \\[dB/m\\]. Dimensões (Nc, Nf).
- `Ti`: Current transformation matrix. Dimensões (Nc, Nc, Nf).

# Notes
The relation between Neper and decibels is `1 dB = log(10) / 20 Np`.
"""
function modos_propagacao(
    Z::Array{Complex{T}, 3},
    Y::Array{Complex{T}, 3},
    complex_frequencies::Vector{Complex{T}},
    unwrap::Bool = true;
    tol::Float64 = 1e-8,
    max_iter::Int = 1000,
) where {T <: Real}

    Nf = length(complex_frequencies)
    Nc = size(Z, 1)
    Np_to_db = log(10) / 20  # conversion factor from Neper to dB
    omega = 2pi * imag.(complex_frequencies)
    w2u0e0 = -omega.^2 * μ₀ * ε₀

    propagation = zeros(Complex{T}, Nc, Nf)
    velocity = zeros(T, Nc, Nf)
    attenuation = zeros(T, Nc, Nf)
    Ti = similar(Z)

    if unwrap
        S_omega = zeros(Complex{T}, Nc, Nc, Nf)
        for i in 1:Nf
            YZ = Y[:, :, i] * Z[:, :, i]
            S_omega[:, :, i] = YZ / (-w2u0e0[i]) - I
        end
        evals, evecs = eig_levenberg_marquardt(S_omega; tol = tol , max_iter = max_iter)
        for i in 1:Nf
            evals[:, i] = (-w2u0e0[i] .* (1 .+ evals[:, i]))
            propagation[:, i] = sqrt.(evals[:, i])
            velocity[:, i] = @. 1e-6 * imag(complex_frequencies[i]) / imag(propagation[:, i])
            attenuation[:, i] = real(propagation[:, i]) * Np_to_db
            Ti[:, :, i] = evecs[:, :, i]
        end
    else
        for i in 1:Nf
            YZ = Y[:, :, i] * Z[:, :, i]
            evals, evecs = eigen(YZ)
            Ti[:, :, i] = evecs
            propagation[:, i] .= sqrt.(evals)
            velocity[:, i] = @. 1e-6 * imag(complex_frequencies[i]) / imag(propagation[:, i])
            attenuation[:, i] = real(propagation[:, i]) * Np_to_db
        end
    end
    return (propagation), (velocity), (attenuation), Ti
end


using NPZ
zc = npzread("test/fixtures/tripolar_Z.npy")
yc = npzread("test/fixtures/tripolar_Y.npy")
gamma_esperado = npzread("test/fixtures/tripolar_propagacao.npy")

# Frequências
nf = 200
freq = exp10.(range(0, 9, nf))
freq_s = freq * 2im * pi

for i in 1:nf
    gamma_esperado[:, i] = sort(gamma_esperado[:, i]; by = cplxpair)
end
for unwrap in [false, true]
    gamma, velocidade, atenuacao, Ti = modos_propagacao(
        zc, yc, freq_s, unwrap, tol = 1e-8, max_iter = 5000
    )
    for i in 1:nf
        gamma[:, i] = sort(gamma[:, i]; by = cplxpair)
    end
    @test gamma ≈ gamma_esperado
end