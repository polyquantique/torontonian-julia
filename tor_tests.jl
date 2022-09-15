#Tests for the Torontonian functions

using DoubleFloats
using IterTools
using LinearAlgebra
using Random
using RandomMatrices
using SpecialFunctions
using Test
include("./torontonian.jl")

"""
    random_covariance(N; hbar=2, pure=false)

Generate random covariance matrix of a Gaussian state.

### Input
- `N`  -- Int: Number of modes
- `hbar` -- Float: (optional, default: `hbar=2`) value of ``\\hbar`` in the commutation relation ``[x,p]=i\\hbar``
- `pure` -- Bool: (optional, default: `pure=false`) if true, generates a matrix corresponding to a pure Gaussian state
### Output

Array: ``2N \\times 2N`` random covariance matrix.
"""
function random_covariance(N; hbar = 2, pure = false)
    U = rand(Haar(2), N)
    O = [real(U) -imag(U); imag(U) real(U)]
    U = rand(Haar(2), N)
    P = [real(U) -imag(U); imag(U) real(U)]
    r = abs.(randn(N))
    sq = diagm(vcat(exp.(-r), exp.(r)))
    S = O * sq * P
    if pure
        return (hbar / 2) .* (S * transpose(S))
    else
        nbar = 2 .* (rand(N)) .+ 1
        V = (hbar / 2) * diagm(vcat(nbar, nbar))
        return S * V * transpose(S)
    end
end

"""
    gen_omats(l, nbar, dtype=ComplexDF64)

Generates the matrix O that enters inside the Torontonian for an l mode system
in which the first mode is prepared in a thermal state with mean photon number nbar
and the rest in vacuum and are later sent into a Fourier interferometer, i.e. one described
by a DFT unitary matrix.

### Input
- `l`    -- Int: number of modes
- `nbar` -- Float: mean photon number of the first mode (the only one not prepared in vacuum)
- `dtype`-- (optional, default: `dtype=ComplexDF64`) type of the output matrix 

### Output

Array: an O matrix whose Torontonian can be calculated analytically.
"""
function gen_omats(l, nbar; dtype = ComplexDF64)
    A = (nbar / (l * (1.0 + nbar))) * ones(dtype, l, l)
    O = [A 0*A; 0*A A]
    return O
end

"""
    analytical_tor(l, nbar)

Return the value of the Torontonian of the O matrices generated by gen_omats.

### Input
- `l`    -- Int: number of modes
- `nbar` -- Float: mean photon number of the first mode (the only one not prepared in vacuum)

### Output

Float: value of the torontonian of gen_omats(l, nbar).
"""
function analytical_tor(l, nbar)
    if isapprox(l, nbar, atol = 1e-14, rtol = 0)
        return 1.0
    end
    beta = -(nbar / (l * (1 + nbar)))
    pref = factorial(l) / beta
    p1 = (pref * l * gamma(1 / beta)) / gamma(1 / beta + l + 2)
    p2 = (pref * beta * gamma(2 + 1 / beta)) / gamma(2 + 1 / beta + l)
    return (p1 + p2) * ((-1)^l)
end

@testset verbose = true "Torontonian tests" begin

    @testset "Torontonian (from definition) vs. analytical result" begin
        for i = 2:1:13
            for k = 1.1:0.37:7.3
                @test isapprox(analytical_tor(i, k), tor(gen_omats(i, k)))
            end
        end
    end

    @testset "Torontonian (recursive) vs. analytical result" begin
        for i = 2:1:13
            for k = 1.1:0.37:7.3
                @test isapprox(analytical_tor(i, k), rec_tor(gen_omats(i, k)))
            end
        end
    end

    @testset "Torontonian (from definition) vs. Torontonian (recursive)" begin
        for i = 2:1:13
            for k = 1.1:0.37:7.3
                @test isapprox(tor(gen_omats(i, k)), rec_tor(gen_omats(i, k)))
            end
        end
    end

    @testset "Torontonian of two-mode squeezed vacuum with mean photon number 1.0" begin
        mean_n = 1.0
        r = asinh(sqrt(mean_n))
        phase = exp(1im * 0.3)
        phasec = conj(phase)
        Omat =
            tanh(r) .*
            Array{ComplexDF64}([0 0 0 phase; 0 0 phase 0; 0 phasec 0 0; phasec 0 0 0])
        @test isapprox(real(tor(Omat)), 1.0)
        @test isapprox(real(rec_tor(Omat)), 1.0)
        @test isapprox(real(tor(Omat)), real(rec_tor(Omat)))
    end

end

@testset "probabilities sum to 1" begin
    for n = 2:1:10
        cov = random_covariance(n)
        prob_df = 0.0
        prob_rc = 0.0
        for p in product([[0 1] for i = 1:n]...)
            pattern = collect(p)
            prob_df += threshold_detection_prob(cov, pattern, recursive = false)
            prob_rc += threshold_detection_prob(cov, pattern, recursive = true)
        end
        @test isapprox(prob_df, 1.0, rtol = 10^(-10))
        @test isapprox(prob_rc, 1.0, rtol = 10^(-10))
    end
end

@testset "probabilities are equal in both methods" begin
    for n = 2:1:10
        cov = random_covariance(n)
        for p in product([[0 1] for i = 1:n]...)
            pattern = collect(p)
            prob_df = threshold_detection_prob(cov, pattern, recursive = false)
            prob_rc = threshold_detection_prob(cov, pattern, recursive = true)
            @test isapprox(prob_df, prob_rc, rtol = 10^(-10))
        end
    end
end
