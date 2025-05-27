function _integrand_real_Lindhard_finT_old(omb, qb, x, bbar)
    num = QEDprobing._nu_minus(omb, qb)
    nup = QEDprobing._nu_plus(omb, qb)
    term1 = log(abs((x - num) / (x + num)))
    term2 = log(abs((x - nup) / (x + nup)))
    return _F(x, bbar) * (term1 - term2) / 2
end

function real_lindhard_nonzero_temperature_old(ombar::T, qbar::T, bbar::T) where {T <: Real}

    num = QEDprobing._nu_minus(ombar, qbar)
    nup = QEDprobing._nu_plus(ombar, qbar)

    res, _ = quadgk(
        x -> _integrand_real_Lindhard_finT_old(ombar, qbar, x, bbar),
        zero(num),
        max(num, zero(num)),
        nup,
        Inf,
    )

    return -res / (qbar)
end
