"""

    _TS_diffCS_pol_spin_sum(cth)

Return unpolarized differential cross section of Thomson scattering as a function of the cosine of the scattering angle.

!!! note "Internal usage"

    This function is only for internal usage. To access the differential cross section of Thomson scattering, use the
    appropriate interface.

"""
function _TS_diffCS_pol_spin_sum(cth)
    return ALPHA_SQUARE / 2 * (1 + cth^2)
    #return ALPHA_SQUARE / 2 * (1 + cth^2)
end

function QEDbase.unsafe_differential_cross_section(
    psp::PhaseSpacePoint{
        QEDprobing.ThomsonScattering,
        PerturbativeQED,
        QEDprobing.ElabPhotonSphSystem,
    },
)
    k_prime = momentum(psp, Outgoing(), Photon())
    cth = getCosTheta(k_prime)
    return QEDprobing._TS_diffCS_pol_spin_sum(cth)
end

@inline function QEDbase._is_in_phasespace(
    psp::PhaseSpacePoint{<:ThomsonScattering,PerturbativeQED},
)
    return isapprox(sum(momenta(psp, Incoming())), sum(momenta(psp, Outgoing())))
end
