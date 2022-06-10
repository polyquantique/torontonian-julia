using DoubleFloats
using Combinatorics
using LinearAlgebra

function q_mat(cov, hbar=2)
    n = size(cov)[1] ÷ 2
    x = (2 / hbar) .* cov[1:n, 1:n]
    xp = (2 / hbar) .* cov[1:n, n + 1:end]
    p = (2 / hbar) .* cov[n + 1:end, n + 1:end]
    ad = (x + p + im .* (xp - transpose(xp)) - 2 .* I) ./ 4
    a = (x - p + im .* (xp + transpose(xp))) ./ 4
    return [ad conj(a); a conj(ad)] + I
end

function tor(O)
    O = Array{ComplexDF64}(O)
    n = size(O)[1] ÷ 2
    total = 0.0
    for set in powerset(1:n)
        pm_coeff = (-1) ^ (length(set))
        kept_rows = vcat(set, set .+ n)
        btt = sqrt(det(I - O[kept_rows, kept_rows]))
        total += pm_coeff / btt
    end
    return ((-1) ^ n) * total
end

function threshold_detection_prob(cov, det_pattern, hbar=2)
    cov = Array{ComplexDF64}(cov)
    n = size(cov)[1] ÷ 2
    qcov = q_mat(cov, hbar)
    O = I - inv(qcov)
    rpatt = vcat(det_pattern, det_pattern)
    rows = findall(x -> x != 0, rpatt)
    rO = O[rows, rows]
    return tor(rO) / sqrt(det(qcov))
end