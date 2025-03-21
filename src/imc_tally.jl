# IMC Tallying Function
# Simon Butson

module Tally
using ..Constants


function tally(inputs, mesh, simvars, particles)
"""Tally end of timestep quantities 
    Parameters:
    inputs: Dict - Dictionary of input parameters
    mesh: Mesh - Mesh object
    simvars: SimVars - Simulation variables object
    particles: Array{Array{Float64,1},1} - Array of particle data
    Returns:
    None
"""

    precision = simvars.precision
    if simvars.geometry == "1D"
        sigma_a = mesh.sigma_a[:,1]
    elseif simvars.geometry == "2D"
        sigma_a = mesh.sigma_a[:,:,1]
    end

    # Radiation energy density
    radenergydens = zeros(precision, mesh.Ncells)
    radenergydens = phys_a * (mesh.temp).^4 # keV/cm^3

    #nrg_inc = (mesh.energydep ./ (mesh.dx .* mesh.dy' * mesh.energyscale )) .- (sigma_a .* mesh.fleck .* radenergydens * phys_c * simvars.dt)
  
    nrg_inc = (mesh.energydep ./ (mesh.dx .* mesh.dy' * mesh.energyscale)) .- mesh.emittedenergy

    #print("The energy increase is ", nrg_inc, "\n")
   
    mesh.matenergydens += nrg_inc

    mesh.temp = mesh.temp .+  nrg_inc ./ mesh.bee

    #print("The updated mesh temperature is ", mesh.temp)

    # Save radiation energy
    mesh.radenergydens = zeros(precision, mesh.Ncells)

    for ii in eachindex(particles)
        if inputs["GEOMETRY"] == "1D"
            rad_index = Int(particles[ii][3])
            mesh.radenergydens[rad_index] += particles[ii][7]/(mesh.dx[rad_index])
        else
            rad_xindex = Int(particles[ii][2])
            rad_yindex = Int(particles[ii][3])
            mesh.radenergydens[rad_xindex,rad_yindex] += particles[ii][8]/(mesh.dx[rad_xindex]*mesh.dx[rad_yindex])
        end
    end
    mesh.radenergydens =  mesh.radenergydens ./ mesh.energyscale

    mesh.totalenergydep += sum(mesh.energydep)

    push!(mesh.temp_saved, copy(mesh.temp)) # Add the current mesh temperature to a saved list

    push!(mesh.matenergy_saved, copy(mesh.matenergydens)) # Add the current material energy to a saved list

    push!(mesh.radenergy_saved, copy(mesh.radenergydens)) # Add the current radiation energy to a saved list


    #total_energyincrease = sum(nrg_inc[:])/ Ncells
    #print("The total energy increase was ", total_energyincrease, "\n")
    #push!(energy_increases, total_energyincrease)

end
end