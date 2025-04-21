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
include("../src/imc_energycheck.jl")
include("../src/imc_utilities.jl")
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
using .EnergyCheck
using .Utilities
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

# @testset "IMC Main Loop Test" begin
#     test_input_path = joinpath(@__DIR__, "test_input.txt")
#     MixedPrecisionIMC.main([test_input_path])
# end

@testset "IMC Update Tests" begin
    test_input_path = joinpath(@__DIR__, "test_input.txt")
    input_data = MixedPrecisionIMC.Input.readInputs(test_input_path)
    mesh = MixedPrecisionIMC.Mesh.mesh_generation(input_data)

    MixedPrecisionIMC.Constants.set_constants(input_data)
    import MixedPrecisionIMC.Constants: phys_c, phys_a, alpha

    simvars = MixedPrecisionIMC.SimVars(0.00, 0.01, 0.0, 0.0, 0.0, 1.0, 0.01, 500, 100000, 1, "FALSE", "REFLECTIVE", Float64, "2D")

    MixedPrecisionIMC.Update.update(mesh, simvars)
    
    @test !iszero(mesh.beta)  # Check that the value of beta is not zero
    @test !iszero(mesh.sigma_a[:,:,1])  # Check that the value of sigma_a is not zero
    @test iszero(mesh.sigma_s[:,:,1])  # Check that the value of sigma_s is zero
    
   
end

@testset "IMC Clean Tests" begin
    # Create a particles list and a single particle to remove
    test_data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, -1.0]
    particles = Vector{Vector{}}()
    push!(particles, test_data)

    @test !isempty(particles)  # Check that particles list is not empty
    MixedPrecisionIMC.Clean.clean(particles)
    @test isempty(particles)  # Check if the cleaned particles list is now empty
end