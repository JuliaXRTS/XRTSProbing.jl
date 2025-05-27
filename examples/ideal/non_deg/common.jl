using Plots
using QuadGK
using LaTeXStrings
using BenchmarkTools
using SpecialFunctions
using Printf
using Plots.PlotMeasures
using LogExpFunctions

import QEDprobing

PLOTDIR = "plots"

function _chemical_potential_normalized(bbar)
    A = 0.25954
    B = 0.072
    b = 0.858

    term1 = 3 / 2 * log(bbar)
    term2 = log(4 / (3 * sqrt(pi)))
    term3 = (A * bbar^(b + 1) + B * bbar^((b + 1) / 2)) / (1 + A * bbar^b)

    return (term1 + term2 + term3) / bbar
end

function general_fermi(x, beta)
    mu = _chemical_potential_normalized(beta)
    denom = exp(beta * (x^2 - mu)) + one(beta)
    return inv(denom)
end

function general_chem_pot_check(beta)
    res, _ = quadgk(x -> general_fermi(x, beta) * x^2, 0, Inf)
    return abs(3 * res - one(beta))
end

function general_integrand(x, nu, beta)
    return x * general_fermi(x, beta) * log(abs(x - nu) / abs(x + nu))
end

function general_integal(nu, beta)
    res, _ = quadgk(x -> general_integrand(x, nu, beta), 0.0, nu, Inf)

    return res
end

function general_integal2(num, nup, beta)
    res, _ = quadgk(
        x -> [general_integrand(x, num, beta), general_integrand(x, nup, beta)],
        0.0,
        num,
        Inf,
    )

    return res
end

function _non_deg_chem_pot(bbar)
    term = 4 * bbar^(3 / 2) / (3 * sqrt(pi))
    return log(term) / bbar
end

function non_deg_fermi1(x, beta)
    mu = _non_deg_chem_pot(beta)
    denom = exp(beta * (x^2 - mu)) + one(beta)
    return inv(denom)
end

function non_deg_chem_pot_check(beta)
    res, _ = quadgk(x -> non_deg_fermi1(x, beta) * x^2, 0, Inf)
    return abs(3 * res - one(res))
end

function non_deg_integal(nu, beta)
    mu = _chemical_potential_normalized(beta)
    #mu = _non_deg_chem_pot(beta)

    prefac = -exp(mu * beta) / beta

    return prefac * dawson(nu * sqrt(beta)) * sqrt(pi)
end

function _real_lindhard_non_deg_old(ombar, qbar, bbar)
    mubar = _non_deg_chem_pot(bbar)
    num = QEDprobing._nu_minus(ombar, qbar)
    nup = QEDprobing._nu_plus(ombar, qbar)

    prefac = sqrt(pi) * exp(mubar * bbar) / (2 * qbar * bbar)
    return prefac * (dawson(num * sqrt(bbar)) - dawson(nup * sqrt(bbar)))
end

function _real_lindhard_non_deg(ombar, qbar, bbar)
    num = QEDprobing._nu_minus(ombar, qbar)
    nup = QEDprobing._nu_plus(ombar, qbar)

    prefac = inv(2 * qbar)
    return -prefac * (non_deg_integal(num, bbar) - non_deg_integal(nup, bbar))
end

function _real_lindhard_new(ombar, qbar, bbar)
    num = QEDprobing._nu_minus(ombar, qbar)
    nup = QEDprobing._nu_plus(ombar, qbar)

    prefac = inv(2 * qbar)

    res = eps()
    try
        res = -prefac * (general_integal(num, bbar) - general_integal(nup, bbar))
    catch
        @warn "problem with $ombar,$qbar,$bbar"
        res = eps()
    end
    return res

end
