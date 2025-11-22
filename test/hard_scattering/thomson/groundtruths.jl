function _groundtruth_TS_diffCS(
        in_pol::P1,
        out_pol::P2,
        cth,
        phi,
    ) where {P1 <: AbstractPolarization, P2 <: AbstractPolarization}
    if abs(cth) == one(cth)
        return _groundtruth_TS_diffCS_forward(in_pol, out_pol)
    end

    return _groundtruth_TS_diffCS_not_forward(in_pol, out_pol, cth, phi)
end

# forward scattering
_groundtruth_TS_diffCS_forward(
    ::P1,
    ::P2,
) where {P1 <: AbstractPolarization, P2 <: AbstractPolarization} = 0.0
_groundtruth_TS_diffCS_forward(::P, ::P) where {P <: AbstractPolarization} =
    ElectronicStructureModels.ALPHA_SQUARE


# not forward scattering

# unpolarized

function _groundtruth_TS_diffCS_not_forward(::AllPol, ::AllPol, cth, phi)
    return ElectronicStructureModels.ALPHA_SQUARE / 2 * (1 + cth^2)
end

# polarized

function _groundtruth_TS_diffCS_not_forward(::PolX, ::PolX, cth, phi)
    pol_fac = (cth * cos(phi))^2
    return ElectronicStructureModels.ALPHA_SQUARE * pol_fac
end

function _groundtruth_TS_diffCS_not_forward(::PolX, ::PolY, cth, phi)
    pol_fac = (sin(phi))^2
    return ElectronicStructureModels.ALPHA_SQUARE * pol_fac
end

function _groundtruth_TS_diffCS_not_forward(::PolY, ::PolX, cth, phi)
    pol_fac = (cth * sin(phi))^2
    return ElectronicStructureModels.ALPHA_SQUARE * pol_fac
end

function _groundtruth_TS_diffCS_not_forward(::PolY, ::PolY, cth, phi)
    pol_fac = (cos(phi))^2
    return ElectronicStructureModels.ALPHA_SQUARE * pol_fac
end
