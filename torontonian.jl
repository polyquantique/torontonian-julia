#Julia implementation of the Torontonian function and threshold detector probabilities (without displacement)

using Combinatorics
using DoubleFloats
using LinearAlgebra

"""
    q_mat(cov; hbar=2)

Compute the Husimi covariance matrix of a Gussian state.

### Input
- `cov`  -- Array: ``2N \\times 2N`` ``xxpp``-Wigner covariance matrix of the Gaussian state
- `hbar` -- Float: (optional, default: `hbar=2`) value of ``\\hbar`` in the commutation relation ``[x,p]=i\\hbar``

### Output

Array: complex Husimi covariance matrix of the Gaussian state.
"""
function q_mat(cov; hbar = 2)
    n = size(cov)[1] ÷ 2
    x = (2 / hbar) .* cov[1:n, 1:n]
    xp = (2 / hbar) .* cov[1:n, n+1:end]
    p = (2 / hbar) .* cov[n+1:end, n+1:end]
    ad = (x + p + im .* (xp - transpose(xp)) - 2 .* I) ./ 4
    a = (x - p + im .* (xp + transpose(xp))) ./ 4
    return [ad conj(a); a conj(ad)] + I
end

"""
    quad_cholesky(L, Z, idx, mat)

Returns the Cholesky factorization of a matrix using sub-matrix of prior
Cholesky based on the new matrix and lower right quadrant.

Algorithm from:
https://arxiv.org/pdf/2109.04528.pdf

### Input
- `L`   -- Array: previous Cholesky factorization
- `Z`   -- Array: new sub-matrix indices
- `idx` -- Int: index of starting row/column of lower right quadrant
- `mat` -- Array: new matrix

### Output

Array: the Cholesky factorization of matrix `mat`.
"""
function quad_cholesky(L, Z, idx, mat)
    L = Array{ComplexDF64}(L)
    mat = Array{ComplexDF64}(mat)
    Ls = L[Z, Z]
    lmat = size(mat)[1]
    for i = idx:1:lmat
        for j = idx:1:i-1
            z = 0.0
            for k = 1:1:j-1
                z += Ls[i, k] * conj(Ls[j, k])
            end
            Ls[i, j] = (mat[i, j] - z) / Ls[j, j]
        end
        z = 0.0
        for k = 1:1:i-1
            z += Ls[i, k] * conj(Ls[i, k])
        end
        Ls[i, i] = real(sqrt(mat[i, i] - z))
    end
    return Ls
end

"""
    recursiveTor(L, modes, A, n)

Returns the recursive Torontonian sub-computation of a matrix.

Algorithm from:
https://arxiv.org/pdf/2109.04528.pdf

### Input
- `L`     -- Array: current Cholesky factorization
- `modes` -- Array: optical modes
- `A`     -- Array: a square, Hermitian array of even dimensions
- `n`     -- Int: size of the original matrix 

### Output

Complex: the recursive Torontonian sub-computation of matrix `A`.
"""
function recursiveTor(L, modes, A, n)
    L = Array{ComplexDF64}(L)
    A = Array{ComplexDF64}(A)
    tot = 0.0
    if length(modes) == 0
        start = 1
    else
        start = modes[end] + 1
    end
    for i = start:1:n
        nextModes = vcat(modes, i)
        nm = size(A)[1] ÷ 2
        idx = 2 * (i - length(modes)) - 1
        Z = vcat(collect(1:1:(idx-1)), collect((idx+2):1:(2*nm)))
        Az = A[Z, Z]
        Ls = quad_cholesky(L, Z, idx, I - Az)
        det = prod(diag(Ls)) .^ 2
        tot += ((-1)^length(nextModes)) / sqrt(det) + recursiveTor(Ls, nextModes, Az, n)
    end
    return tot
end

"""
    tor(O)

Returns the Torontonian of a matrix (using directly the definition).

Definition from:
https://arxiv.org/abs/1807.01639

### Input
- `O` -- Array: a square, Hermitian array of even dimensions.

### Output

Complex: the Torontonian of matrix `O`.
"""
function tor(O)
    O = Array{ComplexDF64}(O)
    n = size(O)[1] ÷ 2
    total = 0.0
    for set in powerset(1:n)
        pm_coeff = (-1)^(length(set))
        kept_rows = vcat(set, set .+ n)
        btt = sqrt(det(I - O[kept_rows, kept_rows]))
        total += pm_coeff / btt
    end
    return ((-1)^n) * total
end

"""
    rec_tor(A)

Returns the Torontonian of a matrix (using the recursive algorithm).

Algorithm from:
https://arxiv.org/pdf/2109.04528.pdf

### Input
- `A` -- Array: a square, Hermitian array of even dimensions.

### Output

Complex: the Torontonian of matrix `A`.
"""
function rec_tor(A)
    A = Array{ComplexDF64}(A)
    n = size(A)[1] ÷ 2
    Z = Array{Int64}(undef, 2 * n)
    Z[1:2:end] = collect(1:1:n)
    Z[2:2:end] = collect((n+1):1:(2*n))
    A = A[Z, Z]
    L = cholesky(I - A).L
    det = prod(diag(L)) .^ 2
    return 1 / sqrt(det) + recursiveTor(L, Array{Int64}(undef, 0), A, n)
end

"""
    threshold_detection_prob(cov, det_pattern; hbar=2, recursive=true)

Returns threshold detection probabilities for Gaussian states without displacement.

### Input
- `cov`         -- Array: 2D xp Wigner covariance matrix
- `det_pattern` -- Array: 1D array of {0,1} that describes the threshold detection outcome
- `hbar`        -- Float: (optional, default: `hbar=2`) value of ``\\hbar`` in the 
                          commutation relation ``[x,p]=i\\hbar``
- `recursive`   -- Bool: (optional, default: `recursive=true`) indicates whether to use the 
                        recursive algorithm for the computation of the Torontonian or not

### Output

Float: probability of the detection pattern.
"""
function threshold_detection_prob(cov, det_pattern; hbar = 2, recursive = true)
    cov = Array{ComplexDF64}(cov)
    n = size(cov)[1] ÷ 2
    qcov = q_mat(cov, hbar = hbar)
    O = I - inv(qcov)
    rpatt = vcat(det_pattern, det_pattern)
    rows = findall(x -> x != 0, rpatt)
    rO = O[rows, rows]
    if recursive
        return rec_tor(Hermitian(rO)) / sqrt(det(qcov))
    else
        return tor(rO) / sqrt(det(qcov))
    end
end
