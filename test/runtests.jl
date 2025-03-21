# Mixed Precision IMC Unit Tests
# Simon Butson

# MixedPrecisionIMC Modules
include("../src/MixedPrecisionIMC.jl")
include("../src/imc_input.jl")
include("../src/constants.jl")
include("../src/imc_mesh.jl")
include("../src/imc_update.jl")
include("../src/imc_sourcing.jl")
include("../src/imc_transport.jl")
include("../src/imc_clean.jl")
include("../src/imc_tally.jl")
include("../src/imc_output.jl")

using MixedPrecisionIMC
using Test
using Random
using .Input
using .Constants
using .Mesh
using .Update
using .Sourcing
using .Transport
using .Clean
using .Tally


@testset "Sample Test Set" begin
    @test 1 + 1 == 2
    @test "hello" == "hello"
end

@testset "IMC Input Tests" begin
    test_input_path = joinpath(@__DIR__, "test_input.txt")
    input_data = MixedPrecisionIMC.Input.readInputs(test_input_path)

    @test !isempty(input_data)  # Check if the input data is not empty
end

@testset "IMC Clean Tests" begin
    # Assuming there is a function `clean_data` in the Clean module
    test_data = [1, 2, 3, 4, 5]
    cleaned_data = MixedPrecisionIMC.Clean.clean_data(test_data)

    @test !isempty(cleaned_data)  # Check if the cleaned data is not empty
    @test length(cleaned_data) == length(test_data)  # Check if the length of cleaned data is the same as the input data
    @test all(x -> x in test_data, cleaned_data)  # Check if all elements in cleaned data are from the input data
end