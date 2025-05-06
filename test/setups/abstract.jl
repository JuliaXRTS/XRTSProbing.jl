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

# setup with default implementations
struct TestSetupDefault <: AbstractTestSetup end

# setup with custom _assert_valid_input
struct TestSetupCustomAssertValidInput <: AbstractTestSetup end
QEDprobing._assert_valid_input(stp::TestSetupCustomAssertValidInput, x) =
    _groundtruth_valid_input_assert(x)

# setup with custom post processing
struct TestSetupCustomPostProcessing <: AbstractTestSetup end
QEDprobing._post_processing(::TestSetupCustomPostProcessing, x, y) =
    _groundtruth_post_processing(x, y)

# setup with custom input validation and post processing
struct TestSetupCustom <: AbstractTestSetup end
QEDprobing._assert_valid_input(stp::TestSetupCustom, x) = _groundtruth_valid_input_assert(x)
QEDprobing._post_processing(::TestSetupCustom, x, y) = _groundtruth_post_processing(x, y)

# setup which fail on computation with default implementations
struct TestSetupFAIL <: AbstractComputationSetup end

# setup which fail on computation with custom input validation, where the
# invalid input will be caught before the computation.
struct TestSetupCustomValidationFAIL <: AbstractComputationSetup end
QEDprobing._assert_valid_input(stp::TestSetupCustomValidationFAIL, x) =
    _groundtruth_valid_input_assert(x)

# setup which fail on computation with custom post processing
struct TestSetupCustomPostProcessingFAIL <: AbstractComputationSetup end
QEDprobing._post_processing(::TestSetupCustomPostProcessingFAIL, x, y) =
    _groundtruth_post_processing(x, y)


@testset "general computation setup interface" begin
    @testset "interface fail" begin
        rnd_input = rand(RNG)

        @test_throws MethodError QEDprobing._compute(TestSetupFAIL(), rnd_input)
        @test_throws MethodError compute(TestSetupFAIL(), rnd_input)

        @test_throws MethodError QEDprobing._compute(
            TestSetupCustomValidationFAIL(),
            rnd_input,
        )
        @test_throws MethodError compute(TestSetupCustomValidationFAIL(), rnd_input)
        # invalid input should be caught without throwing a MethodError
        @test_throws TestException compute(
            TestSetupCustomValidationFAIL(),
            _transform_to_invalid(rnd_input),
        )


        @test_throws MethodError QEDprobing._compute(
            TestSetupCustomPostProcessingFAIL(),
            rnd_input,
        )
        @test_throws MethodError compute(TestSetupCustomPostProcessingFAIL(), rnd_input)
    end

    @testset "default interface" begin
        stp = TestSetupDefault()

        rnd_input = rand(RNG)
        rnd_output = rand(RNG)
        @test QEDprobing._post_processing(stp, rnd_input, rnd_output) == rnd_output
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

    @testset "custom input validation" begin
        stp = TestSetupCustomAssertValidInput()
        rnd_input = rand(RNG)
        @test QEDprobing._assert_valid_input(stp, rnd_input) == nothing
        @test_throws TestException QEDprobing._assert_valid_input(
            stp,
            _transform_to_invalid(rnd_input),
        )
        @test_throws TestException compute(stp, _transform_to_invalid(rnd_input))

    end

    @testset "custom post processing" begin
        stp = TestSetupCustomPostProcessing()
        rnd_input = rand(RNG)
        rnd_output = rand(RNG)
        @test isapprox(
            QEDprobing._post_processing(stp, rnd_input, rnd_output),
            _groundtruth_post_processing(rnd_input, rnd_output),
        )
        @test isapprox(
            compute(stp, rnd_input),
            _groundtruth_post_processing(rnd_input, _groundtruth_compute(rnd_input)),
        )
    end

    @testset "custom input validation and post processing" begin
        stp = TestSetupCustom()
        rnd_input = rand(RNG)
        rnd_output = rand(RNG)

        @test_throws TestException() compute(stp, _transform_to_invalid(rnd_input))
        @test isapprox(
            QEDprobing._post_processing(stp, rnd_input, rnd_output),
            _groundtruth_post_processing(rnd_input, rnd_output),
        )
        @test isapprox(
            compute(stp, rnd_input),
            _groundtruth_post_processing(rnd_input, _groundtruth_compute(rnd_input)),
        )
    end
end


# process setup

const N_INCOMING = rand(RNG, 2:8)
const N_OUTGOING = rand(RNG, 2:8)

const test_process = MockProcess(RNG, N_INCOMING, N_OUTGOING)
const test_model = MockModel()
const test_psl = MockOutPhaseSpaceLayout(MockMomentum)

struct TestProcessSetup <: AbstractProcessSetup end
QEDbase.process(::TestProcessSetup) = test_process
QEDbase.model(::TestProcessSetup) = test_model
QEDbase.phase_space_layout(::TestProcessSetup) = test_psl

struct TestProcessSetupFAIL <: AbstractProcessSetup end

@testset "process setup interface" begin
    @testset "interface fail" begin
        rnd_input = rand(RNG)
        @test_throws MethodError process(TestProcessSetupFAIL())
        @test_throws MethodError model(TestProcessSetupFAIL())
        @test_throws MethodError QEDprobing._compute(TestProcessSetupFAIL(), rnd_input)
        @test_throws MethodError compute(TestProcessSetupFAIL(), rnd_input)
    end

    @testset "hard interface" begin
        stp = TestProcessSetup()

        @test process(stp) == test_process
        @test model(stp) == test_model
        @test phase_space_layout(stp) == test_psl
    end

    @testset "delegated functions" begin
        stp = TestProcessSetup()
        @test number_incoming_particles(stp) == N_INCOMING
        @test number_outgoing_particles(stp) == N_OUTGOING
    end

end
