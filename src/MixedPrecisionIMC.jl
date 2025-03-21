# Main Function for Mixed Precision IMC Software Package
# Simon Butson

module MixedPrecisionIMC

using Random
#using Plots

include("imc_input.jl")
include("constants.jl")
include("imc_mesh.jl")
include("imc_update.jl")
include("imc_sourcing.jl")
include("imc_transport.jl")
include("imc_clean.jl")
include("imc_tally.jl")
include("imc_energycheck.jl")
include("imc_output.jl")

using .Input
using .Constants
using .Mesh
using .Update
using .Sourcing
using .Transport
using .Clean
using .Tally
using .EnergyCheck

# Main Control loop

mutable struct SimVars
    t
    dt
    dt0
    k
    dtmax
    t_end
    timesteps
    n_input::Int
    n_max::Int
    BC
    precision::DataType
    geometry
end

mutable struct RWVars
    aVals
    prVals
    ptVals
end

function main()
    """ Main function to run the Mixed Precision IMC simulation
    Parameters:
    None
    Returns:
    None
    """

    input_file = ARGS[1] # Pass the input file as an argument from the command line

    #input_file = raw"src\inputs\SuOlson.txt"

    inputs = Input.readInputs(input_file)

    Constants.set_constants(inputs)

    mesh = Mesh.mesh_generation(inputs)

    print(inputs, "\n")

    
    Random.seed!(parse(Int, inputs["SEED"]))

    precision = inputs["PRECISION"]
    timestepping = uppercase(string(inputs["TIMESTEPPING"]))


    if timestepping == "CONSTANT"
        dt = parse(precision,inputs["DT"])
        t_end = parse(precision, inputs["ENDTIME"])
        dt0 = precision(0.0)
        k = precision(0.0)
        dtmax = precision(0.0)
    elseif timestepping == "RAMP"
        dt0 = parse(precision,inputs["DT0"])
        k = parse(precision,inputs["K"])
        dtmax = parse(precision,inputs["DTMAX"]) 
        t_end = parse(precision, inputs["ENDTIME"])
        dt = dt0
    end

    n_input = parse(Int, inputs["NINPUT"])
    n_max = parse(Int, inputs["NMAX"])
    geometry = inputs["GEOMETRY"]

    t = 0.0
    timesteps = Vector{precision}()

    particles = Vector{Vector{}}()

    if geometry == "1D"

        LBC = uppercase(string(inputs["LEFTBC"]))
        RBC = uppercase(string(inputs["RIGHTBC"]))
        BC = (LBC, RBC)

        if inputs["RANDOMWALK"] == "TRUE"
            # Create lookup table from probabilities
            aVals = precision.(LinRange(0,10,1000))
            prVals = zeros(precision, 1000)
            ptVals = zeros(precision, 1000)
            prVals, ptVals = Transport.randomwalk_table(aVals, prVals, ptVals)
            rwvars = RWVars(aVals, prVals, ptVals)
        end

        simvars = SimVars(t, dt, dt0, k, dtmax, t_end, timesteps, n_input, n_max, BC, precision, geometry)
        while simvars.t <= simvars.t_end
            print("Time: ", simvars.t, "\n")
            Update.update(mesh, simvars)
            Sourcing.sourcing(mesh, simvars, particles)
            if inputs["RANDOMWALK"] == "TRUE"
                Transport.MC_RW(mesh, simvars, rwvars, particles)
            else
            Transport.MC(mesh, simvars, particles)
            end
            Clean.clean(particles)
            Tally.tally(inputs, mesh, simvars, particles)
            EnergyCheck.energychecker(inputs, mesh, particles)
            Output.plotting(inputs, mesh, simvars)
            timestep(timestepping, simvars)
        end
    elseif geometry == "2D"

        LBC = uppercase(string(inputs["LEFTBC"]))
        RBC = uppercase(string(inputs["RIGHTBC"]))
        TBC = uppercase(string(inputs["TOPBC"]))
        BBC = uppercase(string(inputs["BOTTOMBC"]))
        BC = (LBC, RBC, TBC, BBC)

        simvars = SimVars(t, dt, dt0, k, dtmax, t_end, timesteps, n_input, n_max, BC, precision, geometry)
        while simvars.t <= simvars.t_end
            print("Time: ", simvars.t, "\n")
            Update.update(mesh, simvars)
            Sourcing.sourcing(mesh, simvars, particles)
            Transport.MC2D(mesh, simvars, particles)
            Clean.clean(particles)
            Tally.tally(inputs, mesh, simvars, particles)
            EnergyCheck.energychecker(inputs, mesh, particles)
            Output.plotting(inputs, mesh, simvars)
            timestep(timestepping, simvars)
        end
    end
   
end

function timestep(timestepping, simvars)
""" Function to update the time and time step for the simulation
    Parameters:
    simvars.t: Float - Current time
    simvars.dt: Float - Current time step
    simvars.dtmax: Float - Maximum time step
    simvars.k: Float - Time step multiplier
    simvars.t_end: Float - End time of simulation
    timestepping: String - Time stepping method
    Returns:
    simvars.t: Float - Updated time
"""
    if simvars.t == simvars.t_end
        simvars.t = simvars.t_end + simvars.dt
        return
    end

    if timestepping == "CONSTANT"
        if simvars.t + simvars.dt > simvars.t_end
            simvars.dt = simvars.t_end - simvars.t
            simvars.t = simvars.t_end
        else
            simvars.t = simvars.t + simvars.dt
        end
        elseif timestepping == "RAMP"
        if simvars.dt < simvars.dtmax
            simvars.dt = simvars.dt * simvars.k
            if simvars.dt > simvars.dtmax
                simvars.dt = simvars.dtmax
            end
        end
        if simvars.t + simvars.dt > simvars.t_end
            simvars.dt = simvars.t_end - simvars.t
            simvars.t = simvars.t_end
        else
        simvars.t = simvars.t + simvars.dt
        end
    end
    push!(simvars.timesteps, simvars.t)
    return
end

runtime = @elapsed begin
    main() 
end

print("The simulation took ", runtime, " seconds to run \n")

end

