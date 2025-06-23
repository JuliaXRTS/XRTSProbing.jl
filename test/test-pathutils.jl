using DrWatson

# needed to live in the right project directory
@quickactivate :QEDprobing

DATADIR = datadir()

@testset "datadir" begin
    @testset "raw" begin
        @test datadir_raw() == joinpath(DATADIR, "raw")
        @test datadir_raw("foo") == joinpath(DATADIR, "raw", "foo")
        @test datadir_raw("foo", "bar") == joinpath(DATADIR, "raw", "foo", "bar")

        @test datadir_raw_ext() == joinpath(DATADIR, "raw", "external")
        @test datadir_raw_ext("foo") == joinpath(DATADIR, "raw", "external", "foo")
        @test datadir_raw_ext("foo", "bar") ==
            joinpath(DATADIR, "raw", "external", "foo", "bar")
    end

    @testset "processed" begin
        @test datadir_processed() == joinpath(DATADIR, "processed")
        @test datadir_processed("foo") == joinpath(DATADIR, "processed", "foo")
        @test datadir_processed("foo", "bar") ==
            joinpath(DATADIR, "processed", "foo", "bar")

        @test datadir_processed_ext() == joinpath(DATADIR, "processed", "external")
        @test datadir_processed_ext("foo") ==
            joinpath(DATADIR, "processed", "external", "foo")
        @test datadir_processed_ext("foo", "bar") ==
            joinpath(DATADIR, "processed", "external", "foo", "bar")
    end
end

@testset "readdir" begin
    @testset "raw" begin
        @test readdir_raw_data() == readdir(joinpath(DATADIR, "raw"))
        @test readdir_raw_ext_data() == readdir(joinpath(DATADIR, "raw", "external"))
    end

    @testset "processed" begin
        @test readdir_processed_data() == readdir(joinpath(DATADIR, "processed"))
        @test readdir_processed_ext_data() ==
            readdir(joinpath(DATADIR, "processed", "external"))
    end
end
