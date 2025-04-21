# IMC Energy Conservation Check Script
# Simon Butson


module EnergyCheck

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
        for ii in eachindex(particles)
            push!(radenergy_vector, particles[ii][7]/particles[ii][9])
        end 
    else
        for ii in eachindex(particles)
            push!(radenergy_vector, particles[ii][8]/particles[ii][10])
        end 
    end

    radenergy = sum(radenergy_vector)


    #print("The total energy lost is: ", mesh.lostenergy, " \n")

    energy_error = (mesh.totalenergy - mesh.totalenergydep - radenergy - mesh.lostenergy)/mesh.totalenergy

    print("The energy conservation error is: ", energy_error, " \n")

end

end