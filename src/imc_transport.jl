# IMC Transport Loop
# Simon Butson

module Transport

include("imc_utilities.jl")

using Random
using Base.Iterators
using .Utilities
using ..Constants

function MC(mesh, simvars, particles)
    """ 1D Monte Carlo Transport Function
    Parameters:
    mesh: Mesh - Mesh object
    simvars: SimVars - Simulation variables object
    particles: Array{Array{Float64,1},1} - Array of particle data
    Returns:
    None
    """

    endsteptime = simvars.dt
    precision = simvars.precision
    pairwise = simvars.pairwise
    distancescale = mesh.distancescale
    dx = mesh.dx
    Ncells = mesh.Ncells
    sigma_a = mesh.sigma_a[:,1]
    sigma_s = mesh.sigma_s[:,1]

    if pairwise == "TRUE"
        energydep_vectors = Vector{Vector{precision}}()

        for ii in 1:mesh.Ncells
            for kk in 1:length(mesh.energyscales)
                push!(energydep_vectors, [])
            end
        end   
    end

    mesh.energydep = zeros(precision, (Ncells, length(mesh.energyscales)))
    for particle in eachindex(particles)
    
        currentparticle = particles[particle, :] 
        #print(currentparticle, "\n")
        origin = currentparticle[1][1]
        currenttime = currentparticle[1][2]
        cellindex = Utilities.tointeger(currentparticle[1][3])
        position = currentparticle[1][4] 
        mu = currentparticle[1][5]
        freq = currentparticle[1][6]
        energy = currentparticle[1][7]
        startenergy = currentparticle[1][8]
        energyscale = currentparticle[1][9]
        energyscaleindex = findfirst(isequal(energyscale), mesh.energyscales)

        minenergy = precision(0.01 * startenergy) # Minimum particle energy cut-off

        while true
            simvars.iterations += 1
            # Calculate distances to boundary, collision, and census -> take minimum

            # Boundary distance
            if mu > 0.0
                #dist_b = (mesh_nodes[cellindex+1]-position)/mu
                dist_b = (dx[cellindex] * distancescale - position)/mu
            else
                #dist_b = (mesh_nodes[cellindex]-position)/mu # negatives cancel out
                dist_b = abs(position/mu)
            end
        
            # Collision distance (Scattering)
            #dist_col = abs(randexp()) / sigma_s
            dist_col = (randexp(precision) / (sigma_a[cellindex] * (1 - mesh.fleck[cellindex]) + sigma_s[cellindex])) # Effective scattering includes isotropic term
            # Census distance
            dist_cen = phys_c * (endsteptime - currenttime) * distancescale

            # Actual distance - closest distnace
            dist = min(dist_b, dist_col, dist_cen)

            # Calculate new energy and deposit lost energy
            newenergy = energy * (exp(-sigma_a[cellindex] * mesh.fleck[cellindex]*dist)) # Separate exponent to prevent overflow in argument from large opacities
            #newenergy = energy * exp(-sigma_a * dist)
            if newenergy <= minenergy
                newenergy = precision(0.0)
            end

            # Deposit the particle's energy
            if pairwise == "TRUE"
                push!(energydep_vectors[(energyscaleindex-1)*mesh.Ncells + cellindex], energy - newenergy)
            else
                mesh.energydep[cellindex, energyscaleindex] += energy - newenergy
            end


            # Advance position, time, and energy
            # If energy is zero or domain boundary crossed -> kill particle
            if newenergy == 0.0
                # Flag particle for later destruction
                particles[particle][8] = -1.0
                break
            end

            # Otherwise, advance the position, time and energy
            position += mu * dist
            currenttime += (dist/distancescale) / phys_c 
            energy = newenergy

            # If the event was a boundary-crossing, and the boundary is the
            # domain boundary, then kill the particle
            if dist == dist_b
                if mu > 0
                    if cellindex == Ncells # Right-boundary
                        if simvars.BC[2] == "REFLECT"
                            # Reflecting boundary
                            mu = -mu
                        elseif simvars.BC[2] == "VACUUM"
                            # Vacuum boundary
                            mesh.lostenergy += energy/energyscale
                            particles[particle][8] = -1.0
                            break
                        end
                    end
                    cellindex += 1
                    position = 0 
                end
                if mu < 0
                    if cellindex == 1 # Left-boundary
                        if simvars.BC[1] == "REFLECT"
                            # Reflecting boundary
                            mu = -mu
                        elseif simvars.BC[1] == "VACUUM"
                            # Vacuum boundary
                            mesh.lostenergy += energy/energyscale
                            particles[particle][8] = -1.0
                            break
                        end
                    else
                        cellindex -= 1
                        position = dx[cellindex]*distancescale
                    end
                end
            end
         

            # If collision occured update frequency and direction -> return to cross-section calculation
            if dist == dist_col
            # Collision (i.e. absorption, but treated as pseudo-scattering)
                # freq = phys_invh * mesh_temp[cellindex] * abs(log(rand()))
                #freq = mesh_temp[cellindex] * randexp()
                mu = precision(0.0)
                while mu == 0.0
                    mu = precision(1 - 2*rand(precision))
                end

                end
            # If census event occured, finish history and update particle properties in list
            if dist == dist_cen
                # Finished with this particle
                # Update the particle's properties in the list
                # Starting energy doesn't change
                currenttime = precision(0.0)
                particles[particle][:] = [origin, currenttime, cellindex, position, mu, freq, energy, startenergy, energyscale]
                #global n_census += 1
                break
            end
        end 
    # New particle history
    #print("Particle state at end of time-step ", particles[particle][:], "\n")
    end
    if pairwise == "TRUE"
        for ii in 1:mesh.Ncells
            for kk in 1:length(mesh.energyscales)
                mesh.energydep[ii,kk] = sum(energydep_vectors[(kk-1)*mesh.Ncells + ii])
            end
        end
    end
    #print(mesh.energydep[1,1], "\n")
    print("There were ", simvars.iterations, " total iterations this time-step. ")
    #print("The number of particles that were absorbed is ", absorbed_particles, "\n")
end

function MC_RW(mesh, simvars, rwvars, particles)
    """ 1D Monte Carlo Random Walk Transport Function
    Parameters:
    mesh: Mesh - Mesh object
    simvars: SimVars - Simulation variables object
    rwvars: RWVars - Random Walk variables object
    particles: Array{Array{Float64,1},1} - Array of particle data
    Returns:
    None
    """

    endsteptime = simvars.dt
    dx = mesh.dx
    Ncells = mesh.Ncells
    precision = simvars.precision
    pairwise = simvars.pairwise
    sigma_a = mesh.sigma_a[:,1]
    sigma_s = mesh.sigma_s[:,1]
    aVals = rwvars.aVals
    prVals = rwvars.prVals
    ptVals = rwvars.ptVals

    if pairwise == "TRUE"
        energydep_vectors = Vector{Vector{precision}}()

        for ii in 1:mesh.Ncells
            for kk in 1:length(mesh.energyscales)
                push!(energydep_vectors, [])
            end
        end   
    end

    mesh.energydep = zeros(precision, (Ncells, length(mesh.energyscales)))
    for particle in eachindex(particles)
    
        currentparticle = particles[particle, :] 
        #print(currentparticle, "\n")
        origin = currentparticle[1][1]
        currenttime = currentparticle[1][2]
        cellindex = Utilities.tointeger(currentparticle[1][3])
        position = currentparticle[1][4] 
        mu = currentparticle[1][5]
        freq = currentparticle[1][6]
        energy = currentparticle[1][7]
        startenergy = currentparticle[1][8]
        energyscale = currentparticle[1][9]
        energyscaleindex = findfirst(isequal(energyscale), mesh.energyscales)

        minenergy = precision(0.01 * startenergy) # Minimum particle energy cut-off

        while true
            simvars.iterations += 1
            # Calculate distances to boundary, collision, and census -> take minimum

            # Boundary distance
            if mu > 0.0
                #dist_b = (mesh_nodes[cellindex+1]-position)/mu
                dist_b = (dx[cellindex] - position)/mu
            else
                #dist_b = (mesh_nodes[cellindex]-position)/mu # negatives cancel out
                dist_b = abs(position/mu)
            end
        
            # Collision distance (Scattering)
            #dist_col = abs(randexp()) / sigma_s
            dist_col = abs(randexp()) / (sigma_a[cellindex] * (1 - mesh.fleck[cellindex]) + sigma_s[cellindex]) # Effective scattering includes isotropic term
            # Census distance
            dist_cen = phys_c * (endsteptime - currenttime) 

            # Actual distance - closest distnace
            dist = min(dist_b, dist_col, dist_cen)

            R0 = min(abs(dx[cellindex]-position), abs(position))

            # Check if particle will undergo random walk and initiate if true
            if R0 > 1/mesh.sigma[cellindex] && dist_col < R0
                u = rand(precision)

                D = precision(phys_c/(3*sigma_a[cellindex]*(1-mesh.fleck[cellindex])))

                a = D*endsteptime/(R0*R0)

                Pr = P_r(a, simvars)
                Pt = 1-Pr
                if u < Pt
                    # Solve for actual time t when particle left

                    a_index = bisection(ptVals,u)

                    t_p = precision(aVals[a_index]*(R0*R0)/D)

                    # Updated Absorption Method for Random Walk (not using small fleck factor approximation)
                    newenergy = energy * exp(t_p*phys_c*(1-mesh.fleck[cellindex])*sigma_a[cellindex]/log(1-mesh.fleck[cellindex]))
                    if newenergy <= startenergy
                        newenergy = 0.0
                    end

                    if isnan(newenergy) == true || isnan(energy) == true
                        print(energy, " ", newenergy, " ", mesh.fleck[cellindex], " ", sigma_a[cellindex], " ", dist, " ", t_p , " ", phys_c*(1-mesh.fleck[cellindex])*sigma_a[cellindex]*t_p/log(1-mesh.fleck[cellindex]),"\n")
                    end

                     # Deposit the particle's energy
                    if pairwise == "TRUE"
                        push!(energydep_vectors[(energyscaleindex-1)*mesh.Ncells + cellindex], energy - newenergy)
                    else
                        mesh.energydep[cellindex, energyscaleindex] += energy - newenergy
                    end

                    if newenergy == 0.0
                    # Flag particle for later destruction
                        particles[particle][8] = -1.0
                        break
                    end
                    # Sample new position off radius of sphere
                    omega = precision(2*pi*rand(precision))
                    position += precision(cos(omega) * R0)

                    currenttime += t_p
                    energy = newenergy

                    mu = precision(cos(omega + rand(precision) * pi)) # Sample new direction using cosine distribution about sphere normal
                    continue
                else
                    #Code for if the particle remains in the random walk sphere at census time
                    u_prime = rand(precision)
                    Pr0 = Pr*u_prime

                    a_index = bisection(prVals, Pr0*u_prime)

                    R1 = sqrt(D*endsteptime/aVals[a_index])

                    # Calculate the new energy and the energy deposited (temp storage)
                    newenergy = energy * exp(phys_c*(1-mesh.fleck[cellindex])*sigma_a[cellindex]*endsteptime/log(1-mesh.fleck[cellindex]))

                    if newenergy <= startenergy
                        newenergy = 0.0
                    end
                    # Deposit the particle's energy
                    if pairwise == "TRUE"
                        push!(energydep_vectors[(energyscaleindex-1)*mesh.Ncells + cellindex], energy - newenergy)
                    else
                        mesh.energydep[cellindex, energyscaleindex] += energy - newenergy
                    end

                    if newenergy == 0.0
                        # Flag particle for later destruction
                        particles[particle][8] = -1.0
                        break
                    end
                    position += cos(2*pi*rand(precision)) * R1
                    mu = precision(1 - 2 * rand(precision))
                    currenttime = precision(0.0)
                    energy = newenergy
                    # Finished with this particle
                    # Update the particle's properties in the list
                    # Starting energy comes in here but doesn't change
                    particles[particle][:] = [origin, currenttime, cellindex, position, mu, freq, energy, startenergy, energyscale]
                    break     
                end
            end

            # Calculate new energy and deposit lost energy
            newenergy = energy * exp(-sigma_a[cellindex] * mesh.fleck[cellindex] * dist)
            if newenergy <= minenergy
                newenergy = precision(0.0)
            end

            # Deposit the particle's energy
            if pairwise == "TRUE"
                push!(energydep_vectors[(energyscaleindex-1)*mesh.Ncells + cellindex], energy - newenergy)
            else
                mesh.energydep[cellindex, energyscaleindex] += energy - newenergy
            end
        
            # Advance position, time, and energy
            # If energy is zero or domain boundary crossed -> kill particle
            if newenergy == 0.0
                # Flag particle for later destruction
                particles[particle][8] = -1.0
                break
            end

            # Otherwise, advance the position, time and energy
            position += mu * dist
            currenttime += dist / phys_c 
            energy = newenergy

            # If the event was a boundary-crossing, and the boundary is the
            # domain boundary, then kill the particle
            if dist == dist_b
                if mu > 0
                    if cellindex == Ncells # Right-boundary
                        if simvars.BC[2] == "REFLECT"
                            # Reflecting boundary
                            mu = -mu
                        elseif simvars.BC[2] == "VACUUM"
                            # Vacuum boundary
                            mesh.lostenergy += energy
                            particles[particle][8] = -1.0
                            break
                        end
                    end
                    cellindex += 1
                    position = precision(0) 
                end
                if mu < 0
                    if cellindex == 1 # Left-boundary
                        if simvars.BC[1] == "REFLECT"
                            # Reflecting boundary
                            mu = -mu
                        elseif simvars.BC[1] == "VACUUM"
                            # Vacuum boundary
                            mesh.lostenergy += energy
                            particles[particle][8] = -1.0
                            break
                        end
                    else
                        cellindex -= 1
                        position = dx[cellindex]
                    end
                end
            end

            # If collision occured update frequency and direction -> return to cross-section calculation
            if dist == dist_col
            # Collision (i.e. absorption, but treated as pseudo-scattering)
                
                mu = precision(1 - 2 * rand(precision))
                while mu == 0.0
                    mu = precision(1 - 2 * rand(precision))
                end
            end
            # If census event occured, finish history and update particle properties in list
            if dist == dist_cen
                # Finished with this particle
                # Update the particle's properties in the list
                # Starting energy doesn't change
                currenttime = precision(0.0)
                particles[particle][:] = [origin, currenttime, cellindex, position, mu, freq, energy, startenergy, energyscale]
                #global n_census += 1
                break
            end
        end 
    # New particle history
    #print("Particle state at end of time-step ", particles[particle][:], "\n")
    end
    if pairwise == "TRUE"
        for ii in 1:mesh.Ncells
            for kk in 1:length(mesh.energyscales)
                mesh.energydep[ii,kk] = sum(energydep_vectors[(kk-1)*mesh.Ncells + ii])
            end
        end
    end

    print("There were ", simvars.iterations, " total iterations this time-step. ")
    #print("The number of particles that were absorbed is ", absorbed_particles, "\n")
end



function MC2D(mesh, simvars, particles)
    """ 2D Monte Carlo Transport Function
    Parameters:
    mesh: Mesh - Mesh object
    simvars: SimVars - Simulation variables object
    particles: Array{Array{Float64,1},1} - Array of particle data
    Returns:
    None
    """

    
    precision = simvars.precision
    distancescale = mesh.distancescale
    pairwise = simvars.pairwise
    sigma_a = mesh.sigma_a[:,:,1]
    sigma_s = mesh.sigma_s[:,:,1]
    mesh.energydep = zeros(precision, (mesh.Ncells[1], mesh.Ncells[2], length(mesh.energyscales)))

    if pairwise == "TRUE"
        energydep_vectors = Vector{Vector{precision}}()

        for ii in 1:mesh.Ncells[1]
            for jj in 1:mesh.Ncells[2]
                for kk in 1:length(mesh.energyscales)
                push!(energydep_vectors, [])
                end
            end
        end   
    end

    for particle in eachindex(particles)
        currentparticle = particles[particle, :] 
        currenttime = currentparticle[1][1]
        xindex = Utilities.tointeger(currentparticle[1][2])
        yindex = Utilities.tointeger(currentparticle[1][3])
        xposition = currentparticle[1][4]
        yposition = currentparticle[1][5]
        mu = currentparticle[1][6]
        frq = currentparticle[1][7]
        energy = currentparticle[1][8]
        startenergy = currentparticle[1][9]
        energyscale = currentparticle[1][10]
        energyscaleindex = findfirst(isequal(energyscale), mesh.energyscales)

        #print("Particle number ", particle, "\n")

        minenergy = precision(0.01 * startenergy)   # Minimum particle energy cut-off
        
        while true
            xvec = precision.([cos(mu),sin(mu)]) # Unit direction vector
            # Calculate distances to boundary, collision, and census -> take minimum

            # Boundary distance
            if xvec[1] > 0
                dist_bx = abs((mesh.dx[xindex]*distancescale - xposition)/xvec[1])
            else
                dist_bx = abs(xposition/xvec[1])
            end

            if xvec[2] > 0
                dist_by = abs((mesh.dy[yindex]*distancescale - yposition)/xvec[2])
            else
                dist_by = abs(yposition/xvec[2])
            end

            # Minimum distance to boundary, ignore NaN's that occur when particle orthogonal to one direction
            if isnan(dist_bx) == true
                dist_b = dist_by
            elseif isnan(dist_by) == true
                dist_b = dist_bx
            else
                dist_b = min(dist_bx, dist_by) 
            end
                        

            # Collision distance
            dist_col = randexp(precision)/(sigma_a[xindex,yindex]*(1-mesh.fleck[xindex,yindex]) + sigma_s[xindex,yindex])
            if dist_col < 0
                print("Warning ", dist_col, " ", sigma_a[xindex,yindex], " ", mesh.fleck[xindex,yindex], "\n")
                print("Temperature ", mesh.temp[xindex,yindex], "\n")
                sleep(1)
            end

            # Census distance
            dist_cen = phys_c * (simvars.dt-currenttime) * distancescale
            
            dist = min(dist_b, dist_col, dist_cen)
            if isnan(dist) == true
                print("Warning ", dist_bx, " ", dist_by, " ", xvec, " ", xindex, " ", yindex, "\n")
                print(dist_b, " ", dist_col, " ", dist_cen, " ", dist, " ", dist_bx, " ", dist_by, "\n")
                sleep(1)
            #print(dist_b, " ", dist_col, " ", dist_cen, " ", dist, " ", dist_bx, " ", dist_by, "\n")
            end

            # Calculate new energy and deposit lost energy
            newenergy = energy * exp(-mesh.fleck[xindex,yindex] * sigma_a[xindex,yindex] * dist)
            if newenergy <= minenergy
                newenergy = precision(0.0)
            end

            # Deposit the particle's energy
            if pairwise == "TRUE"
                if energy < 1e-6 && sigma_a[xindex,yindex] < 0.0
                    push!(energydep_vectors[((xindex-1)*mesh.Ncells[2] + yindex) + (energyscaleindex-1)*mesh.Ncells[1]*mesh.Ncells[2]], energy)  
                    newenergy = precision(0.0)
                    # print(energy, " ", newenergy, " ", energy-newenergy, "\n")
                    # sleep(0.1)
                else
                    push!(energydep_vectors[((xindex-1)*mesh.Ncells[2] + yindex) + (energyscaleindex-1)*mesh.Ncells[1]*mesh.Ncells[2]], energy - newenergy)    
                end
            else
                mesh.energydep[xindex,yindex, energyscaleindex] += energy - newenergy
            end

            if isnan(newenergy) == true || isnan(energy) == true
                print(energy, " ", mesh.fleck[xindex,yindex], " ", sigma_a[xindex,yindex], " ", dist, "\n")
            end

            # Advance position, time, and energy
            # If energy is zero or domain boundary crossed -> kill particle
            if iszero(newenergy) == true 
                # Flag particle for later destruction
                particles[particle][8] = -1.0
                break
            end


            # Otherwise, advance the position, time and energy
            xposition += dist * xvec[1]
            yposition += dist * xvec[2]
            currenttime += (dist/distancescale) / phys_c    
            energy = newenergy   

            # Check if particle has left the cell
            if  dist == dist_bx || dist == dist_by
                if dist_bx < dist_by
                    if cos(mu) > 0
                        if xindex == mesh.Ncells[1] # Right-boundary
                            if simvars.BC[2] == "REFLECT"
                                v = [1,0]
                                xreflected = xvec .- v .* (2(v'xvec))
                                mu = atan(xreflected[2],xreflected[1])
                            elseif simvars.BC[2] == "VACUUM"
                                mesh.lostenergy += energy/energyscale
                                particles[particle][8] = -1.0
                                break
                            end
                            continue
                        end
                        xindex += 1
                        xposition = precision(0.0)
                    else
                        if xindex == 1 # Left-boundary
                            if simvars.BC[1] == "REFLECT"
                                v = [1,0]
                                xreflected = xvec .- v .* (2(v'xvec))
                                mu = atan(xreflected[2],xreflected[1])
                            elseif simvars.BC[1] == "VACUUM"
                                mesh.lostenergy += energy/energyscale
                                particles[particle][8] = -1.0
                                break
                            end
                            continue
                        end
                        xindex -= 1
                        xposition = mesh.dx[xindex]*distancescale
                    end
                else
                    if sin(mu) > 0
                        if yindex == mesh.Ncells[2] # Top-boundary
                            if simvars.BC[3] == "REFLECT"
                                v = [0,1]
                                xreflected = xvec .- v .* (2(v'xvec))
                                mu = atan(xreflected[2],xreflected[1])
                            elseif simvars.BC[3] == "VACUUM"
                                mesh.lostenergy += energy/energyscale
                                particles[particle][8] = -1.0
                                break
                            end    
                            continue
                        end
                        yindex += 1
                        yposition = precision(0.0)
                    else
                        if yindex == 1 # Bottom-boundary
                            if simvars.BC[4] == "REFLECT"
                                v = [0,1]
                                xreflected = xvec .- v .* (2(v'xvec))
                                mu = atan(xreflected[2],xreflected[1])
                            elseif simvars.BC[4] == "VACUUM"
                                mesh.lostenergy += energy/energyscale
                                particles[particle][8] = -1.0
                                break
                            end
                            continue
                        end
                        yindex -= 1
                        yposition = mesh.dy[yindex]*distancescale
                    end
                end
                continue
            end

            if dist == dist_col
                # Isotropic Scattering
                mu = precision(2*pi*rand(precision))
                #print(xindex, " ", yindex, " ", xposition, " ", yposition, " ", currenttime, " ", dist, "\n")
            end

            if dist == dist_cen
                currenttime = precision(0)
                particles[particle][:] = [currenttime, xindex, yindex, xposition, yposition, mu, frq, energy, startenergy, energyscale]
                #print("Done with particle \n")
                break
            end

        end
    end
    if pairwise == "TRUE"
        for ii in 1:mesh.Ncells[1]
            for jj in 1:mesh.Ncells[2]
                for kk in 1:length(mesh.energyscales)
                    mesh.energydep[ii,jj,kk] = sum(energydep_vectors[((ii-1)*mesh.Ncells[2] + jj) + (kk-1)*mesh.Ncells[1]*mesh.Ncells[2]])
                end
            end
        end   
    end
end

function P_r(a, simvars)
    """ Probability of a particle not reaching the surface of sphere with radius R_0
        Parameters:
        R_0: Radius of Sphere
        D:  Diffusion Constant {D = c/(3*(1-f)*sigma_R)}
        t: Time 
        a = D*t/(R_0^2)
        simvars: Simulation variables object
    """
    
    Pr = simvars.precision(0.0)
    if a != 0 
        for n in range(1,100)
            Pr += (-1)^(n-1) * exp(-a*(pi*n)^2)*2
        end
    else
        Pr = simvars.precision(1.0)
    end

    return Pr
end

function bisection(array,value)
    """Given an array, and given a value, returns an index j such that value is between array[j]
    and array[j+1]. The array must be monotonic increasing. j=-1 or j=len(array) is returned
    to indicate that value is out of range below and above respectively."""
    n = length(array)
    if (value < array[1])
        return 1
    elseif (value > array[n])
        return n
    end
    jl = 1# Initialize lower
    ju = n# and upper limits.
    while (ju-jl > 1)# If we are not yet done,
        jm=(ju+jl) >> 1# compute a midpoint with a bitshift
        if (value >= array[jm])
            jl=jm# and replace either the lower limit
        else
            ju=jm# or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
        end
    end
    if (value == array[1])# edge cases at bottom
        return 1
    elseif (value == array[n]) # and top
        return n
    else
        return jl
    end
end

function randomwalk_table(aVals, prVals, ptVals, simvars)
    """Create a lookup table for the probability of a particle not reaching the surface of a sphere
    with radius R_0. The table is created for a range of a values, and the corresponding Pr and Pt values
    are stored in the arrays prVals and ptVals respectively."""
  # Create lookup table
  for i in eachindex(aVals)
    prVals[i] = P_r(aVals[i], simvars)
    ptVals[i] = simvars.precision(1-prVals[i])
  end
    return prVals, ptVals  

end

end