#= Functions for modal analysis. =#

using LinearAlgebra
using LeastSquaresOptim


function make_residuals(S_re, S_im, N)
    return function f(x)
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
end


"""Solve the complex symmetric eigenvalue problem for multiple frequencies
preserving the modal ordering.

S_omega may need to be normalized. Example:
```
w2u0e0 = -omega.^2 * mu_0 * epsilon_0
S_omega = [A[i] / (-w2u0e0[i]) - I for i in 1:K]

eigvals, eigvecs = eig_levenberg_marquardt(S_omega, tol=1e-11, max_iter=10000)

# scale back the eigenvalues; eigenvectors remain unchanged
eigvals = [(-w2u0e0 .* (1 .+ eigvals[:, i])) for i in 1:N]
```

# Parameters
- `S_omega`: (N, N, K) array of complex matrices S(ω) for K frequencies
- `tol`: Tolerance for convergence
- `max_iter`: Maximum number of iterations
- `ignore_error`: throw an Error if convergence is not achieved?
- `optimizer`: optimizer and solver from LeastSquaresOptim, see: https://github.com/matthieugomez/LeastSquaresOptim.jl/tree/main?tab=readme-ov-file#choice-of-optimizer--least-square-solver

# Returns
- `eigenvalues`: (N, K) array of complex eigenvalues
- `eigenvectors`: (N, N, K) array of complex eigenvectors (columns)

# References
- A. I. Chrysochos, T. A. Papadopoulos and G. K. Papagiannis, "Robust Calculation
    of Frequency-Dependent Transmission-Line Transformation Matrices Using the
    Levenberg–Marquardt Method," in IEEE Transactions on Power Delivery, vol. 29,
    no. 4, pp. 1621-1629, Aug. 2014, doi: 10.1109/TPWRD.2013.2284504.

- M. L. A. Lourakis and A. A. Argyros, "Is Levenberg-Marquardt the most efficient
    optimization algorithm for implementing bundle adjustment?," Tenth IEEE
    International Conference on Computer Vision (ICCV'05) Volume 1, Beijing,
    China, 2005, pp. 1526-1531 Vol. 2, doi: 10.1109/ICCV.2005.128.
"""
function sorted_eigen(
    S_omega;
    tol = 1e-12,
    max_iter = 1000,
    ignore_error = false,
    optimizer = Dogleg(LeastSquaresOptim.QR()),
)
    _, N, K = size(S_omega)
    eigenvalues = zeros(ComplexF64, N, K)
    eigenvectors = zeros(ComplexF64, N, N, K)

    # First solve the first frequency with standard eigensolver
    S0 = S_omega[:, :, 1]
    F = eigen(S0; permute = false, scale = true, sortby = nothing)
    eigenvalues[:, 1] = F.values
    eigenvectors[:, :, 1] = F.vectors

    for k in 2:K
        # Prepare real-valued formulation
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
            residuals = make_residuals(S_re, S_im, N)

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

            res = optimize(
                residuals,
                x0,
                optimizer;
                x_tol = tol,
                f_tol = tol,
                g_tol = tol,
                iterations = max_iter,
            )
            converged = res.converged
            res_x = res.minimizer

            if !converged
                msg = "Did not converge for S_omega[:, :, $(k)], eigenvalue $(i)"
                if ignore_error
                    @warn msg
                else
                    @error msg
                end
            end

            # Extract solution
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
- `ignore_error`: throw an Error if convergence is not achieved?

# Returns

- `propagation`: Propagation modes `a + jb`, where `a` is \\[Np/m\\] and `b` is \\[rad/s\\]. Dimensões (Nf, Nc).
- `velocity`: Velocity of the propagation modes in \\[m/μs\\]. Dimensões (Nf, Nc).
- `attenuation`: Attenuation of the propagation modes in \\[dB/m\\]. Dimensões (Nf, Nc).
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
    ignore_error = false,
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
        evals, evecs = sorted_eigen(
            S_omega; tol = tol , max_iter = max_iter, ignore_error = ignore_error
        )

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
    return transpose(propagation), transpose(velocity), transpose(attenuation), Ti
end
