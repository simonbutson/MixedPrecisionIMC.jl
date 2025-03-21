# IMC Energy Conservation Check Script
# Simon Butson


module EnergyCheck

function energychecker(inputs, mesh, particles)
    """ Function to check the energy conservation of the IMC simulation
        Parameters:
        inputs: Dict - Dictionary of input parameters
        mesh: Mesh - Mesh object
        particles: Array{Array{Float64,1},1} - Array of particle data
        Returns:
        energy_error: Float - Energy conservation error
    """
    radenergy = 0
    if inputs["GEOMETRY"] == "1D"
        for ii in eachindex(particles)
            radenergy += particles[ii][7]
        end 
    else
        for ii in eachindex(particles)
            radenergy += particles[ii][8]
        end 
    end

    energy_error = (mesh.totalenergy - mesh.totalenergydep - radenergy - mesh.lostenergy)/mesh.totalenergy

    print("The energy conservation error is: ", energy_error, " \n")

end

end