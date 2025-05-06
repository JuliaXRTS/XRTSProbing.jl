@inline _fermi_wave_vector(ne::T) where {T <: Real} = cbrt(3 * pi^2 * ne)
@inline fermi_wave_vector(esys::AbstractElectronSystem) =
    _fermi_wave_vector(electron_density(esys))

@inline _fermi_energy_from_kF(kF::T) where {T <: Real} = kF^2 / 2
@inline fermi_energy(esys::AbstractElectronSystem) =
    _fermi_energy_from_kF(fermi_wave_vector(esys))

@inline beta(esys::AbstractElectronSystem) = inv(temperature(esys))
@inline betabar(esys::AbstractElectronSystem) = fermi_energy(esys) * beta(esys)

function _transform_om_q(elsys::AbstractElectronSystem, om_q::NTuple{2, T}) where {T <: Real}
    kF = fermi_wave_vector(elsys)
    EF = fermi_energy(elsys)

    om, q = om_q
    ombar = om / EF
    qbar = q / kF

    return ombar, qbar
end

function dynamic_response(esys::AbstractElectronSystem, om_q::NTuple{2, T}) where {T <: Real}
    return real_dynamic_response(esys, om_q) + 1im * imag_dynamic_response(esys, om_q)
end


@inline function _dynamic_structure_factor_pos_om(
        esys::AbstractElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    _fac = pi * electron_density(esys)
    imag_rf = imag_dynamic_response(esys, om_q)

    om = @inbounds om_q[1]
    if iszero(temperature(esys))
        return -imag_rf / _fac
    end

    return -imag_rf / ((one(T) - exp(-beta(esys) * om)) * _fac)
end

function dynamic_structure_factor(
        esys::AbstractElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    om, q = @inbounds om_q
    if om < zero(om)
        return exp(om * beta(esys)) * _dynamic_structure_factor_pos_om(esys, (-om, q))
    else
        return _dynamic_structure_factor_pos_om(esys, om_q)
    end
end

function static_response(esys::AbstractElectronSystem, q::T) where {T <: Real}

end

function static_structure_factor(esys::AbstractElectronSystem, q::T) where {T <: Real}

end
