# finite temperature

function _imag_lindhard_large_om(omb, qb, bbar)
    num = QEDprobing._nu_minus(omb, qb)
    prefac = -pi / (4 * qb * bbar)
    mubar = QEDprobing._chemical_potential_normalized(bbar)
    return prefac * exp(-bbar * (num^2 - mubar))
end

function _om_crit_func(qbar, bbar)
    e = 1.0e-8
    mubar = QEDprobing._chemical_potential_normalized(bbar)
    log_arg = 4 * e / pi * qbar * bbar
    r = abs(mubar - inv(bbar) * log(log_arg))
    return 2 * qbar * sqrt(r) + qbar^2
end

function _imag_lindhard_finT_unstable(ombar, qbar, bbar)
    num = QEDprobing._nu_minus(ombar, qbar)
    nup = QEDprobing._nu_plus(ombar, qbar)
    mubar = _chemical_potential_normalized(bbar)

    term1 = log1pexp(bbar * (num^2 - mubar))
    term2 = log1pexp(bbar * (nup^2 - mubar))

    return -pi / (4 * qbar) * (ombar + (term1 - term2) / bbar)
end

imag_lindhard_nonzero_temperature(ombar::T, qbar::T, bbar::T) where {T <: Real} = imag_lindhard_nonzero_temperature(ombar, qbar, bbar, NoApprox())
function imag_lindhard_nonzero_temperature(ombar::T, qbar::T, bbar::T, ::NoApprox) where {T <: Real}
    if bbar < 5.0 * one(T) &&
            zero(qbar) < abs(qbar) <= 2.0e-3 &&
            abs(ombar) >= _om_crit_func(qbar, bbar)
        return _imag_lindhard_large_om(ombar, qbar, bbar)
    else
        return _imag_lindhard_finT_unstable(ombar, qbar, bbar)
    end
end

## real lindhard -> stable implementation


# this is slighly better than LogExpFunctions version
function _f_prime(xi, A, B)
    #u = A * xi^2 - B
    #sig = logistic(u)
    #return -2 * A * xi * sig*(one(sig)-sig)
    u = A * xi^2 - B
    # Compute sigmoid derivative stably:
    if u ≥ 0
        exp_neg_u = exp(-u)
        s_deriv = exp_neg_u / (1 + exp_neg_u)^2
    else
        exp_u = exp(u)
        s_deriv = exp_u / (1 + exp_u)^2
    end
    return -2 * A * xi * s_deriv
end

function transformed_integrand(xi, A, B)
    f_prime_val = _f_prime(xi, A, B)
    part1 = xi * f_prime_val
    prefac = (xi^2 - 1) / 2
    part2 = -f_prime_val * prefac * log(
        abs(
            (xi - 1) / (xi + 1)
        )
    )

    return part1 + part2
end

function transformed_integral_positive(nu, bbar)
    A = bbar * nu^2
    B = _chemical_potential_normalized(bbar) * bbar
    res, _ = quadgk(x -> transformed_integrand(x, A, B), 0.0, 1.0, Inf)

    return nu^2 * res
end

function transformed_integral(nu, bbar)
    if nu < zero(nu)
        return -transformed_integral_positive(-nu, bbar)
    end
    return transformed_integral_positive(nu, bbar)
end

real_lindhard_nonzero_temperature(ombar::T, qbar::T, bbar::T) where {T <: Real} = real_lindhard_nonzero_temperature(ombar, qbar, bbar, NoApprox())
function real_lindhard_nonzero_temperature(ombar::T, qbar::T, bbar::T, ::NoApprox) where {T <: Real}
    num = QEDprobing._nu_minus(ombar, qbar)
    nup = QEDprobing._nu_plus(ombar, qbar)
    return -(transformed_integral(num, bbar) - transformed_integral(nup, bbar)) / (2 * qbar)
end

#function imag_lindhard_nonzero_temperature(ombar::T, qbar::T, bbar::T) where {T <: Real}


# zero temperature approximation

function _imag_lindhardDSF_zeroT_minus(omb::T, qb::T) where {T <: Real}
    num = _nu_minus(omb, qb)

    if abs(num) >= one(T)
        return zero(T)
    end
    return one(omb) - num^2

end

function _imag_lindhardDSF_zeroT_plus(omb::T, qb::T) where {T <: Real}
    nup = _nu_plus(omb, qb)
    if abs(nup) >= one(T)
        return zero(T)
    end
    return one(omb) - nup^2

end

imag_lindhard_zero_temperature(ombar::T, qbar::T) where {T <: Real} = imag_lindhard_zero_temperature(ombar, qbar, NoApprox())
function imag_lindhard_zero_temperature(ombar::T, qbar::T, ::NoApprox) where {T <: Real}
    return -pi / (4 * qbar) *
        (_imag_lindhardDSF_zeroT_minus(ombar, qbar) - _imag_lindhardDSF_zeroT_plus(ombar, qbar))
end

real_lindhard_zero_temperature(ombar::T, qbar::T) where {T <: Real} = real_lindhard_zero_temperature(ombar, qbar, NoApprox())
function real_lindhard_zero_temperature(ombar::T, qbar::T, ::NoApprox) where {T <: Real}
    num = _nu_minus(ombar, qbar)
    nup = _nu_plus(ombar, qbar)

    term1 = (1 - num^2) / 4 * _stable_log_term(num)
    term2 = (1 - nup^2) / 4 * _stable_log_term(nup)

    return -(qbar / 2 - term1 + term2) / qbar
end
