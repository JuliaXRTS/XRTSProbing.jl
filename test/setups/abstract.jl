using QEDprobing
using QEDbase
using QEDbase.Mocks
using QEDcore

using Random
RNG = MersenneTwister(137137)
ATOL = 0.0
RTOL = sqrt(eps())

include("groundtruths.jl")

# setups for which the interface is implemented
abstract type AbstractTestSetup <: AbstractComputationSetup end
QEDprobing._compute(stp::AbstractTestSetup, x) = _groundtruth_compute(x)

# test setup
struct TestSetup <: AbstractTestSetup end

# setup which fail on computation with default implementations
struct TestSetupFAIL <: AbstractComputationSetup end

@testset "general computation setup interface" begin
    @testset "interface fail" begin
        rnd_input = rand(RNG)

        @test_throws MethodError QEDprobing._compute(TestSetupFAIL(), rnd_input)
        @test_throws MethodError compute(TestSetupFAIL(), rnd_input)
    end

    @testset "interface" begin
        stp = TestSetup()

        rnd_input = rand(RNG)
        @test isapprox(
            QEDprobing._compute(stp, rnd_input),
            _groundtruth_compute(rnd_input),
            atol = ATOL,
            rtol = RTOL,
        )
        @test isapprox(
            compute(stp, rnd_input),
            _groundtruth_compute(rnd_input),
            atol = ATOL,
            rtol = RTOL,
        )
    end
end


# process setup

const N_INCOMING = rand(RNG, 2:8)
const N_OUTGOING = rand(RNG, 2:8)

const test_process = MockProcess(RNG, N_INCOMING, N_OUTGOING)
const test_model = MockModel()
const test_psl = MockOutPhaseSpaceLayout(MockMomentum)

struct TestProcessSetup{P, M, PSL} <: AbstractProcessSetup{P, M, PSL}
    proc::P
    model::M
    psl::PSL
end
QEDbase.process(stp::TestProcessSetup) = stp.proc
QEDbase.model(stp::TestProcessSetup) = stp.model
QEDbase.phase_space_layout(stp::TestProcessSetup) = stp.psl

struct TestProcessSetupFAIL{P, M, PSL} <: AbstractProcessSetup{P, M, PSL}
    proc::P
    model::M
    psl::PSL
end

@testset "process setup interface" begin
    @testset "interface fail" begin
        stp = TestProcessSetupFAIL(test_process, test_model, test_psl)
        rnd_input = rand(RNG)
        @test_throws MethodError process(TestProcessSetupFAIL())
        @test_throws MethodError model(TestProcessSetupFAIL())
        @test_throws MethodError QEDprobing._compute(TestProcessSetupFAIL(), rnd_input)
        @test_throws MethodError compute(TestProcessSetupFAIL(), rnd_input)
    end

    @testset "hard interface" begin
        stp = TestProcessSetup(test_process, test_model, test_psl)

        @test process(stp) == test_process
        @test model(stp) == test_model
        @test phase_space_layout(stp) == test_psl
    end

    @testset "delegated functions" begin
        stp = TestProcessSetup(test_process, test_model, test_psl)
        @test number_incoming_particles(stp) == N_INCOMING
        @test number_outgoing_particles(stp) == N_OUTGOING
    end

    # TODO:
    # - implement tests for `_compute`, `compute` for different inputs
    # - implement tests for `_build_psp`
    # - don't forget checking non-valid inputs

end
