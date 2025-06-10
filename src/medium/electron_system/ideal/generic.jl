function imag_dynamic_response(
        elsys::AbstractIdealElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)
    approx = response_approximation(elsys)

    if ombar <= zero(ombar)
        if qbar < zero(qbar)
            return _imag_ideal_dynamic_response(elsys, approx, -ombar, -qbar)
        else
            return -_imag_ideal_dynamic_response(elsys, approx, -ombar, qbar)
        end
    end

    return _imag_ideal_dynamic_response(elsys, approx, ombar, qbar)
end

function real_dynamic_response(
        elsys::AbstractIdealElectronSystem,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)

    # TODO: insert symmetries here as well
    return _real_ideal_dynamic_response(elsys, response_approximation(elsys), ombar, qbar)
end
