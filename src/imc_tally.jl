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
    #radenergydens = zeros(precision, mesh.Ncells) # keV/cm^3

    #nrg_inc = (mesh.energydep ./ (mesh.dx .* mesh.dy' * mesh.energyscale )) .- (sigma_a .* mesh.fleck .* radenergydens * phys_c * simvars.dt)
    
    #print("Emitted energy ", sum(mesh.emittedenergy), "\n")

    nrg_inc = zeros(precision, mesh.Ncells)
    for ii in eachindex(mesh.energyscales)
        if simvars.geometry == "1D"
            #push!(nrg_incvecs, (mesh.energydep[:,ii] ./ (mesh.dx) .- mesh.emittedenergy[:,ii]) / mesh.energyscales[ii])
            nrg_inc += (mesh.energydep[:,ii] ./ (mesh.dx) .- mesh.emittedenergy[:,ii]) / mesh.energyscales[ii]
            mesh.totalenergydep += sum(mesh.energydep[:,ii] / mesh.energyscales[ii])
        elseif simvars.geometry == "2D"
            nrg_inc += (mesh.energydep[:,:,ii] ./ (mesh.dx .* mesh.dy') .- mesh.emittedenergy[:,:,ii]) / mesh.energyscales[ii]
            print("Energy increase in group ", ii, " is ", sum((mesh.energydep[:,:,ii] ./ (mesh.dx .* mesh.dy') - mesh.emittedenergy[:,:,ii]) / mesh.energyscales[ii]), "\n")
            mesh.totalenergydep += sum(mesh.energydep[:,:,ii] / mesh.energyscales[ii])
        end
    end
    push!(mesh.energyincrease_saved, nrg_inc) # Add the current energy increase to a saved list
    #nrg_inc = sum(nrg_incvecs) # Sum over the energy groups

    #nrg_inc = (mesh.energydep ./ (mesh.dx .* mesh.dy') - mesh.emittedenergy)

    #print("Nan's energydep ", sum(isnan.(mesh.energydep)), "\n")
    #print("Nan's nrg_inc ", sum(isnan.(nrg_inc)), "\n")
    #print("Nan's emittedenergy ", sum(isnan.(mesh.emittedenergy)), "\n")
   
    print("Energy increase: ", sum(nrg_inc), "\n")
    mesh.matenergydens += nrg_inc


    mesh.temp = mesh.temp .+  nrg_inc ./ (mesh.bee)
    #mesh.temp = mesh.temp_saved[1] .+ sum(mesh.energyincrease_saved) ./ mesh.bee # Update the temperature based on the energy increase


    #print("The updated mesh temperature is ", mesh.temp)

    # Save radiation energy
    mesh.radenergydens = zeros(precision, mesh.Ncells)


    radenergydens_vectors = Vector{Vector{precision}}()

    # if inputs["GEOMETRY"] == "1D"
    #     for ii in 1:mesh.Ncells
    #         push!(radenergydens_vectors, [])
    #     end   
    #     for jj in eachindex(particles)
    #         rad_index = Int(particles[jj][3])
    #         push!(radenergydens_vectors[rad_index], particles[jj][7]/(mesh.dx[rad_index]*particles[jj][9]))
    #     end
    #     for kk in 1:mesh.Ncells
    #         mesh.radenergydens[kk] = sum(radenergy_vectors[kk])
    #     end
    # elseif inputs["GEOMETRY"] == "2D"
    #     for ii in 1:mesh.Ncells[1]
    #         for jj in 1:mesh.Ncells[2]
    #             push!(radenergydens_vectors, [])
    #         end
    #     end   
    #     for kk in eachindex(particles)
    #         rad_xindex = Int(particles[kk][2])
    #         rad_yindex = Int(particles[kk][3])
    #         push!(radenergydens_vectors[(rad_xindex-1)*mesh.Ncells[2] + rad_yindex], particles[kk][8]/(mesh.dx[rad_xindex]*mesh.dy[rad_yindex]*particles[kk][10]))
    #     end
    #     for ii in 1:mesh.Ncells[1]
    #         for jj in 1:mesh.Ncells[2]
    #             mesh.radenergydens[ii,jj] = sum(radenergydens_vectors[(ii-1)*mesh.Ncells[2] + jj])
    #         end
    #     end   
    # end
 
    for ii in eachindex(particles)
        if inputs["GEOMETRY"] == "1D"
            rad_index = Int(particles[ii][3])
            mesh.radenergydens[rad_index] += particles[ii][7]/(mesh.dx[rad_index]*particles[ii][9])
        else
            rad_xindex = Int(particles[ii][2])
            rad_yindex = Int(particles[ii][3])
            mesh.radenergydens[rad_xindex,rad_yindex] += particles[ii][8]/(mesh.dx[rad_xindex]*mesh.dy[rad_yindex]*particles[ii][10])
        end
    end

    

    push!(mesh.temp_saved, copy(mesh.temp)) # Add the current mesh temperature to a saved list

    push!(mesh.matenergy_saved, copy(mesh.matenergydens)) # Add the current material energy to a saved list

    push!(mesh.radenergy_saved, copy(mesh.radenergydens)) # Add the current radiation energy to a saved list


    #total_energyincrease = sum(nrg_inc[:])/ Ncells
    #print("The total energy increase was ", total_energyincrease, "\n")
    #push!(energy_increases, total_energyincrease)

end
end