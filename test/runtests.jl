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

using Test
using Random
using .Input
using .Mesh
using .Update
using .Sourcing
using .Transport
using .Clean
using .Tally
using ..Constants

@testset "Sample Test Set" begin
    @test 1 + 1 == 2
    @test "hello" == "hello"
end

@testset "IMC Inputs and Mesh Generation Tests" begin
    test_input_path = joinpath(@__DIR__, "test_input.txt")
    input_data = MixedPrecisionIMC.Input.readInputs(test_input_path)

    @test !isempty(input_data)  # Check if the input data is not empty
    @test isa(input_data, Dict)  # Check if the input data is a dictionary
    @test haskey(input_data, "PRECISION")  # Check if the input data has the key "PRECISION"
    
    MixedPrecisionIMC.Constants.set_constants(input_data)
    import MixedPrecisionIMC.Constants: phys_c
    @test isa(MixedPrecisionIMC.Constants.phys_c, Float64)  # Check if phys_c has been set to a Float64 value

    mesh = MixedPrecisionIMC.Mesh.mesh_generation(input_data)
    #@test !isempty(mesh)  # Check if the mesh is not empty
    @test isa(mesh.nodes, Tuple{Vector{Float64}, Vector{Float64}})  # Check if the value of the key "XMESH" is an array
end

# @testset "IMC Clean Tests" begin
#     # Assuming there is a function `clean_data` in the Clean module
#     test_data = [1, 2, 3, 4, 5]
#     cleaned_data = MixedPrecisionIMC.Clean.clean_data(test_data)

#     @test !isempty(cleaned_data)  # Check if the cleaned data is not empty
#     @test length(cleaned_data) == length(test_data)  # Check if the length of cleaned data is the same as the input data
#     @test all(x -> x in test_data, cleaned_data)  # Check if all elements in cleaned data are from the input data
# end