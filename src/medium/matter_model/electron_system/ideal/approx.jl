abstract type AbstractApproximation end

abstract type AbstractResponseApproximation <: AbstractApproximation end

# most naive implementation
struct NoApprox <: AbstractResponseApproximation end

# using integral transforms and identities to make integration more stable
struct TransformedIntegral <: AbstractResponseApproximation end

# using the limit for non-degenerated electron gases
struct NonDegenerated <: AbstractResponseApproximation end


const AbstractZeroTemperatureApproximation = Union{NoApprox}
const AbstractFiniteTemperatureApproximation = Union{NoApprox}
