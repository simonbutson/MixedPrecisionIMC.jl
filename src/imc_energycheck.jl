# IMC Energy Conservation Check Script
# Simon Butson



module EnergyCheck

using ..Constants

function energychecker(inputs, mesh, simvars, particles)
    """ Function to check the energy conservation of the IMC simulation
        Parameters:
        inputs: Dict - Dictionary of input parameters
        mesh: Mesh - Mesh object
        particles: Array{Array{Float64,1},1} - Array of particle data
        Returns:
        energy_error: Float - Energy conservation error
    """
    radenergy_vector = Vector{simvars.precision}()
    if inputs["GEOMETRY"] == "1D"
        # for ii in eachindex(mesh.temp)
        #    radenergy_vector = mesh.radenergydens .* mesh.dx
        # end 
        radenergy_vector = mesh.radenergydens .* mesh.dx
    else
        radenergy_vector = mesh.radenergydens .* mesh.dx .* mesh.dy'
    end

    radenergy = sum(radenergy_vector)

    print("Total energy: ", mesh.totalenergy, " Total energy deposition: ", mesh.totalenergydep, " Radiation energy change: ", radenergy - mesh.radenergyold, " Lost energy: ", mesh.lostenergy, "\n")
    #print("The total energy lost is: ", mesh.lostenergy, " \n")

    energy_error = (mesh.totalenergy - mesh.totalenergydep - (radenergy - mesh.radenergyold) - mesh.lostenergy)/mesh.totalenergy

    mesh.radenergyold = radenergy
    mesh.lostenergy = simvars.precision(0.0)  # Reset the lost energy for the next timestep
    print("The energy conservation error is: ", energy_error, " \n")

    # if abs(energy_error) > 2e-2
    #     print("Warning: Energy conservation error is greater than 2e-2! \n")
    # end

    # Larsen Mercier Time-Step Check
    # T_U = maximum(mesh.temp)
    # max_dt = simvars.precision(0.0)
    # for i in eachindex(mesh.temp)
    #     dt = mesh.bee[i]/(mesh.sigma_a[i] * phys_a * phys_c * ((T_U^4 - mesh.temp[i]^4)/(T_U - mesh.temp[i]) - 4 * alpha * mesh.temp[i]^3)) 
    #     if dt > max_dt
    #         max_dt = dt
    #     end
    # end
    # print("The Larsen Mercier max time-step is: ", max_dt, " \n")


end

end