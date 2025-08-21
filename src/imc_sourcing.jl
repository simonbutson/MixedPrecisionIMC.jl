# IMC Sourcing Function
# Simon Butson

module Sourcing

include("imc_utilities.jl")
using ..Constants
using .Utilities
using Random


    function sourcing(mesh, simvars, particles)
        """
        Function to source particles into the mesh
        Parameters:
        mesh: Mesh - Mesh object    
        simvars: SimVars - Simulation variables object
        particles: Array{Array{Float64,1},1} - Array of particle data
        Returns:    
        None
        """
  
    energyscales = mesh.energyscales
    distancescale = mesh.distancescale
    dt = simvars.dt
    cellmin = simvars.cellmin
    n_input = simvars.n_input
    n_max = simvars.n_max
    T_surface = mesh.temp_surf
    precision = simvars.precision
    Ncells = mesh.Ncells
    e_body = zeros(precision, Ncells)
    e_radsource = zeros(precision, Ncells)
    escale_body = ones(precision, Ncells)
    escale_radsource = ones(precision, Ncells)
    escale_emittedenergy = ones(precision, Ncells)

    if simvars.geometry == "1D"
        sigma_a = mesh.sigma_a[:,1]
    elseif simvars.geometry == "2D"
        sigma_a = mesh.sigma_a[:,:,1]
    end
        # Calculate sourcing energy in each cell

    #print(mesh.temp, "\n")

    e_surfacebottom = zeros(precision, size(mesh.dx))
    e_surfacetop = zeros(precision, size(mesh.dx))
    e_surfaceleft = zeros(precision, size(mesh.dy))
    e_surfaceright = zeros(precision, size(mesh.dy))
    escale_surfacebottom = ones(precision, size(mesh.dx))
    escale_surfacetop = ones(precision, size(mesh.dx))
    escale_surfaceleft = ones(precision, size(mesh.dy))
    escale_surfaceright = ones(precision, size(mesh.dy))

    if simvars.geometry == "1D"

        mesh.emittedenergy = zeros(precision, (Ncells, length(energyscales)))

        (e_surfaceleft, escale_surfaceleft, escaleindex) = Utilities.sorter([phys_a, phys_c, T_surface[1], T_surface[1], T_surface[1], T_surface[1], dt, 0.25], energyscales, precision)
        (e_surfaceright, escale_surfaceright, escaleindex) = Utilities.sorter([phys_a, phys_c, T_surface[2], T_surface[2], T_surface[2], T_surface[2], dt, 0.25], energyscales, precision)
        #e_surfaceleft = precision(energyscale* phys_a * phys_c * (T_surface[1]^4) * dt /4)
        #e_surfaceright = precision(energyscale* phys_a * phys_c * (T_surface[2]^4) * dt /4)
        e_surface = (e_surfaceleft/escale_surfaceleft) + (e_surfaceright/escale_surfaceright)

        for i in eachindex(mesh.temp)
            (e_body[i], escale_body[i], escaleindex) = Utilities.sorter([mesh.fleck[i], sigma_a[i], phys_a, phys_c, mesh.temp[i], mesh.temp[i], mesh.temp[i], mesh.temp[i], mesh.dx[i], dt, distancescale], energyscales, precision)
            (e_radsource[i], escale_radsource[i], escaleindex) = Utilities.sorter([mesh.radsource[i], mesh.dx[i], dt], energyscales, precision)
            (cell_emittedenergy, escale_emittedenergy[i], escaleindex) = Utilities.sorter([mesh.fleck[i], sigma_a[i], phys_a, phys_c, mesh.temp[i], mesh.temp[i], mesh.temp[i], mesh.temp[i], dt, distancescale], energyscales, precision)
            mesh.emittedenergy[i,escaleindex] = cell_emittedenergy
        end
        #print("Body energy ", e_body, "\n")
        #print("Surface source energy ", e_surface, "\n")
        #print("Emitted energy density ", mesh.emittedenergy./escale_emittedenergy, " Total emitted energy density ", sum(mesh.emittedenergy./escale_emittedenergy), "\n")
        print("Total intial time-step energy ", sum(mesh.emittedenergy./escale_emittedenergy) + sum(mesh.radenergydens), "\n")

    elseif simvars.geometry == "2D"

        #e_surfacebottom = precision.(energyscale .* mesh.dx * phys_a * phys_c .* (T_surface[1].^4) * dt /4)
        #e_surfacetop = precision.(energyscale .* mesh.dx * phys_a * phys_c .* (T_surface[2].^4) * dt /4)
        #e_surfaceleft = precision.(energyscale .* mesh.dy * phys_a * phys_c .* (T_surface[3].^4) * dt /4)
        #e_surfaceright = precision.(energyscale .*mesh.dy * phys_a * phys_c .* (T_surface[4].^4) * dt /4)

        mesh.emittedenergy = zeros(precision, (Ncells[1], Ncells[2], length(energyscales)))        

        for i in eachindex(mesh.dx)
            (e_surfacebottom[i], escale_surfacebottom[i], escaleindex) = Utilities.sorter([phys_a, phys_c, T_surface[1][i], T_surface[1][i], T_surface[1][i], T_surface[1][i], mesh.dx[i], dt, 0.25], energyscales, precision)
            (e_surfacetop[i], escale_surfacetop[i], escaleindex) = Utilities.sorter([phys_a, phys_c, T_surface[2][i], T_surface[2][i], T_surface[2][i], T_surface[2][i], mesh.dx[i], dt, 0.25], energyscales, precision)
        end
        for j in eachindex(mesh.dy)
            (e_surfaceleft[j], escale_surfaceleft[j], escaleindex) = Utilities.sorter([phys_a, phys_c, T_surface[3][j], T_surface[3][j], T_surface[3][j], T_surface[3][j], mesh.dy[j], dt, 0.25], energyscales, precision)
            (e_surfaceright[j], escale_surfaceright[j], escaleindex) = Utilities.sorter([phys_a, phys_c, T_surface[4][j], T_surface[4][j], T_surface[4][j], T_surface[4][j], mesh.dy[j], dt, 0.25], energyscales, precision)
        end

        e_surface = sum(e_surfacebottom./escale_surfacebottom) + sum(e_surfacetop./escale_surfacetop) + sum(e_surfaceleft./escale_surfaceleft) + sum(e_surfaceright./escale_surfaceright)
    
        for idx in CartesianIndices(mesh.temp)
            xindex, yindex = Tuple(idx)
            #e_body[xindex,yindex] = precision(mesh.fleck[xindex,yindex] * sigma_a[xindex,yindex] * phys_a * phys_c * (mesh.temp[xindex,yindex]^2) * (mesh.dx[xindex] * mesh.dy[yindex]) * dt * (mesh.temp[xindex,yindex]^2) * distancescale * energyscale) 
            #e_radsource[xindex,yindex] = precision(energyscale * mesh.radsource[xindex,yindex] * (mesh.dx[xindex] * mesh.dy[yindex]) * dt)
            #mesh.emittedenergy[xindex,yindex] = precision(mesh.fleck[xindex,yindex] * sigma_a[xindex,yindex] * phys_a * phys_c * dt * (mesh.temp[xindex,yindex]^4) * distancescale * energyscale)
            (e_body[xindex,yindex], escale_body[xindex,yindex], escaleindex) = Utilities.sorter([mesh.fleck[xindex, yindex], sigma_a[xindex, yindex], phys_a, phys_c, mesh.temp[xindex, yindex], mesh.temp[xindex, yindex], mesh.temp[xindex, yindex], mesh.temp[xindex, yindex], mesh.dx[xindex], mesh.dy[yindex], dt, distancescale], energyscales, precision)
            (e_radsource[xindex,yindex], escale_radsource[xindex,yindex], escaleindex) = Utilities.sorter([mesh.radsource[xindex, yindex], mesh.dx[xindex] , mesh.dy[yindex], dt], energyscales, precision)
            (cell_emittedenergy, escale_emittedenergy[xindex,yindex], escaleindex) = Utilities.sorter([mesh.fleck[xindex, yindex], sigma_a[xindex, yindex], phys_a, phys_c, mesh.temp[xindex, yindex], mesh.temp[xindex, yindex], mesh.temp[xindex, yindex], mesh.temp[xindex, yindex], dt, distancescale], energyscales, precision)
            mesh.emittedenergy[xindex,yindex,escaleindex] = cell_emittedenergy
            #mesh.emittedenergy[xindex,yindex] = e_body[xindex,yindex]/(mesh.dx[xindex] * mesh.dy[yindex])
        end
        #print("Body energy ", sum(e_body), "\n")
        #print("Source energy ", sum(e_radsource), "\n")
        #print("Emitted energy ", sum(mesh.emittedenergy), "\n")
    end 


    #e_body = precision.(mesh.fleck .* sigma_a * phys_a * phys_c .* (mesh.temp.^2) .* (mesh.dx .* mesh.dy') * dt .* (mesh.temp.^2) * distancescale * energyscale) 

    #print("Body energy ", sum(e_body), "\n")

    #e_radsource = precision.(energyscale * mesh.radsource .* (mesh.dx .* mesh.dy') * dt)
    #print("Source energy ", sum(e_radsource), "\n")

    mesh.totalenergy = sum(e_body./escale_body) + sum(e_radsource./escale_radsource) + e_surface   # Total Energy

    #mesh.emittedenergy .= precision.((mesh.fleck .* sigma_a * phys_a * phys_c * dt .* (mesh.temp.^4)) * distancescale * energyscale)
    #print("Emitted energy ", sum(mesh.emittedenergy), "\n")

    # Calculate the number of particles to be sourced

    n_body = zeros(precision, size(e_body))

    n_radsource = zeros(precision, size(e_radsource))

    n_source = n_input
    n_census = length(particles)
    if n_input + n_census > n_max
        n_source = max(cellmin, n_max - n_census - length(Ncells) - 1)
    end

    for cellindex in eachindex(e_body)
        n_body[cellindex] = Utilities.tointeger(max(round((e_body[cellindex]/escale_body[cellindex])*n_source/mesh.totalenergy), cellmin))
    end

    for cellindex in eachindex(e_radsource)
        if e_radsource[cellindex] > 0
            n_radsource[cellindex] = Utilities.tointeger(max(round((e_radsource[cellindex]/escale_radsource[cellindex])*n_source/mesh.totalenergy), cellmin))
        end
    end


    if simvars.geometry == "1D"
        n_surfleft = Utilities.tointeger(precision(0))
        if e_surfaceleft > 0
            n_surfleft = Utilities.tointeger(round(precision, (e_surfaceleft/escale_surfaceleft) * n_source / mesh.totalenergy))
        end
        n_surfright = Utilities.tointeger(precision(0))
        if e_surfaceright > 0
            n_surfright = Utilities.tointeger(round(precision, (e_surfaceright/escale_surfaceright) *n_source / mesh.totalenergy))
        end

        # Surface-source particles
        # Left Surface
        for _ in 1:n_surfleft
            origin = Utilities.tointeger(precision(1))
            xpos = precision(0.01*mesh.dx[1]*distancescale)
            nrg = precision(e_surfaceleft / n_surfleft)
            startnrg = nrg
            mu = precision(sqrt(rand(precision)))
            while mu == 0.0
                mu = precision(sqrt(rand(precision)))
            end
            spawntime = dt * rand(precision)
            #frq = T_surface[1] * sample_planck()
            frq = precision(1.0)
            push!(particles, [origin, spawntime, origin, xpos, mu, frq, nrg, startnrg, escale_surfaceleft])
        end
        # Right Surface
        for _ in 1:n_surfright
            origin = Ncells
            xpos = precision(0.99*mesh.dx[Ncells]*distancescale)
            nrg = precision(e_surfaceright / n_surfright)
            startnrg = nrg
            mu = precision(-sqrt(rand(precision)))
            while mu == 0.0
                mu = precision(-sqrt(rand(precision)))
            end
            spawntime = dt * rand(precision)
            #frq = T_surface[2] * sample_planck()
            frq = precision(1.0)
            push!(particles, [origin, spawntime, origin, xpos, mu, frq, nrg, startnrg, escale_surfaceright])
        end

        # Body-emitted particles
        particle_energysum = precision(0)
        for cellindex in eachindex(e_body)
            if n_body[cellindex] <= 0
                continue
            end
            nrg = e_body[cellindex] / precision(n_body[cellindex])
            startnrg = nrg
            #print("Particle energy ", nrg, "\n")
            for _ in 1:n_body[cellindex]
                particle_energysum += nrg
                origin = cellindex
                xpos = mesh.dx[cellindex] * rand(precision) * distancescale
                mu = precision(1 - 2*rand(precision))
                while mu == 0.0
                    mu = precision(1 - 2*rand(precision))
                end
                spawntime = dt * rand(precision)
                # frq = mesh.temp[cellindex] * sample_planck()
                frq = precision(1.0)
                push!(particles, [origin, spawntime, cellindex, xpos, mu, frq, nrg, startnrg, escale_body[cellindex]])
            end
            #print("Emitted particle energy: ", particle_energysum, "\n")
        end

        # Radiation-source particles
        for cellindex in eachindex(e_radsource)
            if n_radsource[cellindex] <= 0
                continue
            end
            nrg = e_radsource[cellindex] / precision(n_radsource[cellindex])
            startnrg = nrg

            for _ in 1:n_radsource[cellindex]
                origin = cellindex
                xpos = mesh.dx[cellindex] * rand(precision) * distancescale
                mu = precision(1 - 2*rand(precision))
                while mu == 0.0
                    mu = precision(1 - 2*rand(precision))
                end
                spawntime = dt * rand(precision)
                # frq = mesh.temp[cellindex] * sample_planck()
                frq = precision(1.0)
                push!(particles, [origin, spawntime, cellindex, xpos, mu, frq, nrg, startnrg, escale_radsource[cellindex]])
            end
        end

    elseif simvars.geometry == "2D"
         # Surface-source particles
         n_surfbottom = zeros(precision, size(e_surfacebottom))
         for cellindex in eachindex(e_surfacebottom)
            if e_surfacebottom[cellindex] > 0
                n_surfbottom[cellindex] = Utilities.tointeger(max(round((e_surfacebottom[cellindex]/escale_surfacebottom[cellindex])*n_source/mesh.totalenergy),cellmin))
            end
        end
         n_surftop = zeros(precision, size(e_surfacetop))
         for cellindex in eachindex(e_surfacetop)
            if e_surfacetop[cellindex] > 0
                n_surftop[cellindex] = Utilities.tointeger(max(round((e_surfacetop[cellindex]/escale_surfacetop[cellindex])*n_source/mesh.totalenergy),cellmin))
            end
        end
         n_surfleft = zeros(precision, size(e_surfaceleft))
         for cellindex in eachindex(e_surfaceleft)
            if e_surfaceleft[cellindex] > 0
                n_surfleft[cellindex] = Utilities.tointeger(max(round((e_surfaceleft[cellindex]/escale_surfaceleft[cellindex])*n_source/mesh.totalenergy),cellmin))
            end
        end
         n_surfright = zeros(precision, size(e_surfaceright))
         for cellindex in eachindex(e_surfaceright)
            if e_surfaceright[cellindex] > 0
                n_surfright[cellindex] = Utilities.tointeger(max(round((e_surfaceright[cellindex]/escale_surfaceright[cellindex])*n_source/mesh.totalenergy),cellmin))
            end
        end
         # Bottom Surface
         for i in eachindex(e_surfacebottom)
            for _ in 1:n_surfbottom[i]
            spawntime = dt*rand(precision)
            xindex = i
            yindex = 1    
            xpos = mesh.dx[i] * rand(precision) * distancescale
            ypos = precision(0.001*mesh.dy[1]*distancescale)
            mu = precision(pi*rand(precision))    
            frq = precision(1.0)
            nrg = e_surfacebottom[i] / n_surfbottom[i]
            startnrg = nrg
            push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg, escale_surfacebottom[i]]) 
            end
        end
        # Top Surface
        for i in eachindex(e_surfacetop)
            for _ in 1:n_surftop[i]
            spawntime = dt*rand(precision)
            xindex = i
            yindex = Ncells[2]    
            xpos = mesh.dx[i] * rand(precision) * distancescale
            ypos = precision(0.999*mesh.dy[Ncells[2]]*distancescale)
            mu = precision(-pi*rand(precision))    
            frq = precision(1.0)
            nrg = e_surfacetop[i] / n_surftop[i]
            startnrg = nrg
            push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg, escale_surfacetop[i]]) 
            end
        end
        # Left Surface
        for j in eachindex(e_surfaceleft)
            for _ in 1:n_surfleft[j]
            spawntime = dt*rand(precision)
            xindex = 1
            yindex = j    
            xpos = precision(0.001*mesh.dx[1]*distancescale)
            ypos = mesh.dx[j] * rand(precision) * distancescale
            mu = precision(pi*(0.5-rand(precision)))   
            frq = precision(1.0) 
            nrg = e_surfaceleft[j] / n_surfleft[j]
            startnrg = nrg
            push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg, escale_surfaceleft[j]]) 
            end
        end
        # Right Surface
        for j in eachindex(e_surfaceright)
            for _ in 1:n_surfright[j]
            spawntime = dt*rand(precision)
            xindex = Ncells[1]
            yindex = j    
            xpos = precision(0.999*mesh.dx[Ncells[1]]*distancescale)
            ypos = mesh.dx[j] * rand(precision) * distancescale
            mu = precision(pi*(0.5+rand(precision)))   
            frq = precision(1.0) 
            nrg = e_surfaceright[j] / n_surfright[j]
            startnrg = nrg
            push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg, escale_surfaceright[j]]) 
            end
        end

        # Body-emitted particles
        for idx in CartesianIndices(e_body)
            xindex, yindex = Tuple(idx)
            if n_body[xindex, yindex] <= 0
                continue
            end
            nrg = e_body[xindex, yindex] / precision(n_body[xindex, yindex])
            startnrg = nrg
            # if nrg <= floatmin(precision)/100
            #     mesh.emittedenergy[xindex, yindex, :] .= precision(0.0)
            #     continue
            # end
            for _ in 1:n_body[xindex, yindex]
                xpos = mesh.dx[xindex] * rand(precision) * distancescale
                ypos = mesh.dy[yindex] * rand(precision) * distancescale
                mu = precision(2*pi*rand(precision))
                spawntime = dt * rand(precision)
                # frq = mesh.temp[cellindex] * sample_planck()
                frq = precision(1.0)
                push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg, escale_body[xindex, yindex]])
            end
        end

          # Radiation-source particles
          for idx in CartesianIndices(e_radsource)
            xindex, yindex = Tuple(idx)
            if n_radsource[xindex, yindex] <= 0
                continue
            end
            nrg = e_radsource[xindex, yindex] / precision(n_radsource[xindex, yindex])
            startnrg = nrg

            for _ in 1:n_body[xindex, yindex]
                xpos = mesh.dx[xindex] * rand(precision) * distancescale
                ypos = mesh.dy[yindex] * rand(precision) * distancescale
                mu = precision(2*pi*rand(precision))
                spawntime = dt * rand(precision)
                # frq = mesh.temp[cellindex] * sample_planck()
                frq = precision(1.0)
                push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg, escale_radsource[xindex, yindex]])
            end
        end
    end

    print("The number of particles after sourcing is ", length(particles), "\n")
    end

    function sample_planck(simvars)
        """Samples from the Planck spectrum
            Parameters:
            simvars: SimVars - Simulation variables object
            Returns:
            Float - Sampled frequency
        """
        # Method from Fleck and Cummings Paper
        precision = simvars.precision
        n = precision(1.0)
        rn1 = rand(precision)
        nsum = precision(1.0)
       

        while true
            if rn1 <= 90.0 * nsum / pi^4
                rn1 = rand(precision) # Not sure if this should be resampled
                rn2 = rand(precision)
                rn3 = rand(precision)
                rn4 = rand(precision)
                sample_freq = precision(-1.0 * log(rn1 * rn2 * rn3 * rn4) / n)
                return sample_freq
            else
            n += precision(1.0)
            nsum += precision(1.0 / n^4)
            end
        end
    end

end