# TODO:
# - insert approximation type
# - treat ZeroTemperature mode as approximation (ZeroTemperatureApproximation)


struct IdealElectronSystem{TM, T, A <: AbstractResponseApproximation} <: AbstractProperElectronSystem
    electron_density::T
    temperature::T
    approx::A

    function IdealElectronSystem{TM}(
            electron_density::T1,
            temperature::T2,
            approx::A
        ) where {
            T <: Real,
            T1 <: Union{T, Quantity{T}},
            T2 <: Union{T, Quantity{T}},
            TM <: AbstractTemperatureMode,
            A <: AbstractResponseApproximation,
        }

        ne_internal = _internalize_density(electron_density)
        temp_internal = _internalize_temperature(temperature)

        return new{TM, T, A}(ne_internal, temp_internal)
    end
end

# TODO: think about making this a default constructor
function IdealElectronSystem{ZeroTemperature}(
        electron_density::T1,
        approx::A
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
        A <: AbstractZeroTemperatureApproximation,
    }
    return IdealElectronSystem{ZeroTemperature}(electron_density, zero(T), approx)
end

function IdealElectronSystem{ZeroTemperature}(
        electron_density::T1,
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
    }
    return IdealElectronSystem{ZeroTemperature}(electron_density, zero(T), NoApprox())
end

function IdealElectronSystem(
        electron_density::T1,
        temperature::T2,
        approx::A
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
        T2 <: Union{T, Quantity{T}},
        A <: AbstractFiniteTemperatureApproximation,
    }
    return IdealElectronSystem{FiniteTemperature}(electron_density, temperature, approx)
end

function IdealElectronSystem{FiniteTemperature}(
        electron_density::T1,
        temperature::T2,
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
        T2 <: Union{T, Quantity{T}},
    }
    return IdealElectronSystem{FiniteTemperature}(electron_density, temperature, NoApprox())
end

function IdealElectronSystem(
        electron_density::T1,
        temperature::T2,
    ) where {
        T <: Real,
        T1 <: Union{T, Quantity{T}},
        T2 <: Union{T, Quantity{T}},
    }
    return IdealElectronSystem{FiniteTemperature}(electron_density, temperature, NoApprox())
end

@inline temperature(elsys::IdealElectronSystem{ZeroTemperature, T}) where {T <: Real} = zero(T)
@inline temperature(elsys::IdealElectronSystem{FiniteTemperature, T}) where {T <: Real} =
    elsys.temperature
@inline electron_density(elsys::IdealElectronSystem) = elsys.electron_density
@inline response_approximation(elsys::IdealElectronSystem) = elsys.approx


function imag_dynamic_response(
        elsys::IdealElectronSystem{ZeroTemperature},
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)
    return imag_lindhard_zero_temperature(ombar, qbar, response_approximation(elsys))
end

function real_dynamic_response(
        elsys::IdealElectronSystem{ZeroTemperature},
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)
    return real_lindhard_zero_temperature(ombar, qbar, response_approximation(elsys))
end

function imag_dynamic_response(
        elsys::IdealElectronSystem{FiniteTemperature},
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)
    approx = response_approximation(elsys)

    if ombar <= zero(ombar)
        if qbar < zero(qbar)
            return imag_lindhard_nonzero_temperature(-ombar, -qbar, betabar(elsys), approx)
        else
            return -imag_lindhard_nonzero_temperature(-ombar, qbar, betabar(elsys), approx)
        end
    end

    return imag_lindhard_nonzero_temperature(ombar, qbar, betabar(elsys), approx)
end

function real_dynamic_response(
        elsys::IdealElectronSystem{FiniteTemperature},
        om_q::NTuple{2, T},
    ) where {T <: Real}
    ombar, qbar = _transform_om_q(elsys, om_q)
    approx = response_approximation(elsys)
    return real_lindhard_nonzero_temperature(ombar, qbar, betabar(elsys), approx)
end
