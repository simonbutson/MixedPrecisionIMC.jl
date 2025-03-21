# IMC Sourcing Function
# Simon Butson

module Sourcing
using ..Constants
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
    energyscale = mesh.energyscale
    distancescale = mesh.distancescale
    dt = simvars.dt
    n_input = simvars.n_input
    n_max = simvars.n_max
    T_surface = mesh.temp_surf
    precision = simvars.precision
    Ncells = mesh.Ncells
    mesh.emittedenergy = zeros(precision, Ncells)
    if simvars.geometry == "1D"
        sigma_a = mesh.sigma_a[:,1]
    elseif simvars.geometry == "2D"
        sigma_a = mesh.sigma_a[:,:,1]
    end
        # Calculate sourcing energy in each cell

    if simvars.geometry == "1D"
        e_surfaceleft = precision(energyscale* phys_a * phys_c * (T_surface[1]^4) * dt /4)
        e_surfaceright = precision(energyscale* phys_a * phys_c * (T_surface[2]^4) * dt /4)
        e_surface = e_surfaceleft + e_surfaceright
    elseif simvars.geometry == "2D"
        e_surfacebottom = precision.(energyscale .* mesh.dx * phys_a * phys_c .* (T_surface[1].^4) * dt /4)
        e_surfacetop = precision.(energyscale .* mesh.dx * phys_a * phys_c .* (T_surface[2].^4) * dt /4)
        e_surfaceleft = precision.(energyscale .* mesh.dy * phys_a * phys_c .* (T_surface[3].^4) * dt /4)
        e_surfaceright = precision.(energyscale .*mesh.dy * phys_a * phys_c .* (T_surface[4].^4) * dt /4)

        e_surface = sum(e_surfacebottom) + sum(e_surfacetop) + sum(e_surfaceleft) + sum(e_surfaceright)
    end

    e_body = precision.(energyscale * mesh.fleck .* sigma_a * phys_a * phys_c .* (mesh.dx .* mesh.dy') * dt .* (mesh.temp.^4) * distancescale) 

    print("Body energy ", e_body, "\n")

    e_radsource = precision.(energyscale * mesh.radsource .* (mesh.dx .* mesh.dy') * dt)

    e_total = sum(e_body) + sum(e_radsource) + e_surface   # Total Energy

    mesh.emittedenergy .= precision.((mesh.fleck .* sigma_a * phys_a * phys_c * dt .* (mesh.temp.^4)) * distancescale) 

    mesh.totalenergy += e_total

    # Calculate the number of particles to be sourced

    n_body = zeros(precision, size(e_body))

    n_radsource = zeros(precision, size(e_radsource))

    n_source = n_input
    n_census = length(particles)
    if n_input + n_census > n_max
        n_source = n_max - n_census - length(Ncells) - 1
    end

    for cellindex in eachindex(e_body)
        n_body[cellindex] = tointeger(max(round(e_body[cellindex]*n_source/e_total), precision(1)))
    end

    for cellindex in eachindex(e_radsource)
        if e_radsource[cellindex] > 0
            n_radsource[cellindex] = tointeger(max(round(e_radsource[cellindex]*n_source/e_total), precision(1)))
        end
    end


    if simvars.geometry == "1D"
        n_surfleft = tointeger(precision(0))
        if e_surfaceleft > 0
            n_surfleft = tointeger(round(precision, n_source * e_surfaceleft / e_total))
        end
        n_surfright = tointeger(precision(0))
        if e_surfaceright > 0
            n_surfright = tointeger(round(precision, n_source * e_surfaceright / e_total))            
        end

        # Surface-source particles
        # Left Surface
        for _ in 1:n_surfleft
            origin = tointeger(precision(1))
            xpos = precision(0.001*mesh.dx[1]*distancescale)
            nrg = precision(e_surfaceleft / n_surfleft)
            startnrg = nrg
            mu = precision(sqrt(rand(precision)))
            while mu == 0.0
                mu = precision(sqrt(rand(precision)))
            end
            spawntime = dt * rand(precision)
            #frq = T_surface[1] * sample_planck()
            frq = precision(1.0)
            push!(particles, [origin, spawntime, origin, xpos, mu, frq, nrg, startnrg])
        end
        # Right Surface
        for _ in 1:n_surfright
            origin = Ncells
            xpos = precision(0.999*mesh.dx[Ncells]*distancescale)
            nrg = precision(e_surfaceright / n_surfright)
            startnrg = nrg
            mu = precision(-sqrt(rand(precision)))
            while mu == 0.0
                mu = precision(-sqrt(rand(precision)))
            end
            spawntime = dt * rand(precision)
            #frq = T_surface[2] * sample_planck()
            frq = precision(1.0)
            push!(particles, [origin, spawntime, origin, xpos, mu, frq, nrg, startnrg])
        end

        # Body-emitted particles
        for cellindex in eachindex(e_body)
            if n_body[cellindex] <= 0
                continue
            end
            nrg = e_body[cellindex] / precision(n_body[cellindex])
            startnrg = nrg

            for _ in 1:n_body[cellindex]
                origin = cellindex
                xpos = mesh.dx[cellindex] * rand(precision) * distancescale
                mu = precision(1 - 2*rand(precision))
                while mu == 0.0
                    mu = precision(1 - 2*rand(precision))
                end
                spawntime = dt * rand(precision)
                # frq = mesh.temp[cellindex] * sample_planck()
                frq = precision(1.0)
                push!(particles, [origin, spawntime, cellindex, xpos, mu, frq, nrg, startnrg])
            end
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
                push!(particles, [origin, spawntime, cellindex, xpos, mu, frq, nrg, startnrg])
            end
        end

    elseif simvars.geometry == "2D"
         # Surface-source particles
         n_surfbottom = zeros(precision, size(e_surfacebottom))
         for cellindex in eachindex(e_surfacebottom)
            if e_surfacebottom[cellindex] > 0
                n_surfbottom[cellindex] = tointeger(max(round(e_surfacebottom[cellindex]*n_source/e_total),1))
            end
        end
         n_surftop = zeros(precision, size(e_surfacetop))
         for cellindex in eachindex(e_surfacetop)
            if e_surfacetop[cellindex] > 0
                n_surftop[cellindex] = tointeger(max(round(e_surfacetop[cellindex]*n_source/e_total),1))
            end
        end
         n_surfleft = zeros(precision, size(e_surfaceleft))
         for cellindex in eachindex(e_surfaceleft)
            if e_surfaceleft[cellindex] > 0
                n_surfleft[cellindex] = tointeger(max(round(e_surfaceleft[cellindex]*n_source/e_total),1))
            end
        end
         n_surfright = zeros(precision, size(e_surfaceright))
         for cellindex in eachindex(e_surfaceright)
            if e_surfaceright[cellindex] > 0
                n_surfright[cellindex] = tointeger(max(round(e_surfaceright[cellindex]*n_source/e_total),1))
            end
        end
         # Bottom Surface
         for i in eachindex(e_surfacebottom)
            for _ in 1:n_surfbottom[i]
            spawntime = dt*rand()
            xindex = i
            yindex = 1    
            xpos = mesh.dx[i] * rand(precision) * distancescale
            ypos = precision(0.001*mesh.dy[1]*distancescale)
            mu = precision(pi*rand(precision))    
            frq = precision(1.0)
            nrg = e_surfacebottom[i] / n_surfbottom[i]
            startnrg = nrg
            push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg]) 
            end
        end
        # Top Surface
        for i in eachindex(e_surfacetop)
            for _ in 1:n_surftop[i]
            spawntime = dt*rand()
            xindex = i
            yindex = Ncells[2]    
            xpos = mesh.dx[i] * rand(precision) * distancescale
            ypos = precision(0.999*mesh.dy[Ncells[2]]*distancescale)
            mu = precision(-pi*rand(precision))    
            frq = precision(1.0)
            nrg = e_surfacetop[i] / n_surftop[i]
            startnrg = nrg
            push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg]) 
            end
        end
        # Left Surface
        for j in eachindex(e_surfaceleft)
            for _ in 1:n_surfleft[j]
            spawntime = dt*rand()
            xindex = 1
            yindex = j    
            xpos = precision(0.001*mesh.dx[1]*distancescale)
            ypos = mesh.dx[j] * rand(precision) * distancescale
            mu = precision(pi*(0.5-rand(precision)))   
            frq = precision(1.0) 
            nrg = e_surfaceleft[j] / n_surfleft[j]
            startnrg = nrg
            push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg]) 
            end
        end
        # Right Surface
        for j in eachindex(e_surfaceright)
            for _ in 1:n_surfright[j]
            spawntime = dt*rand()
            xindex = Ncells[1]
            yindex = j    
            xpos = precision(0.999*mesh.dx[Ncells[1]]*distancescale)
            ypos = mesh.dx[j] * rand(precision) * distancescale
            mu = precision(pi*(0.5+rand(precision)))   
            frq = precision(1.0) 
            nrg = e_surfaceright[j] / n_surfright[j]
            startnrg = nrg
            push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg]) 
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

            for _ in 1:n_body[xindex, yindex]
                xpos = mesh.dx[xindex] * rand(precision) * distancescale
                ypos = mesh.dy[yindex] * rand(precision) * distancescale
                mu = precision(2*pi*rand(precision))
                spawntime = dt * rand(precision)
                # frq = mesh.temp[cellindex] * sample_planck()
                frq = precision(1.0)
                push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg])
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
                push!(particles, [spawntime, xindex, yindex, xpos, ypos, mu, frq, nrg, startnrg])
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

    function tointeger(x::AbstractFloat)
        """ Function to convert a float to an integer of the same precision
            Parameters:
            x: AbstractFloat - Float to be converted
            Returns:
            Int - Integer value of x with the same precision
        """

        if typeof(x) == Float16
            return Int16(x)
        elseif typeof(x) == Float32
            return Int32(x)
        else
            return Int64(x)
        end
    end
end