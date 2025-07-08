# IMC Mesh Quantity Generation Function
# Simon Butson

module Mesh

include("imc_utilities.jl")
using .Utilities
#using ..Constants

    mutable struct MeshStruct
    temp::Array # Temperature (keV)
    temp_surf # Surface Temperature (keV) Left Boundary
    fleck::Array # Fleck factor
    beta::Array # Beta factor
    sigma::Array # Total opacity
    sigma_a::Array # Absorption opacity
    sigma_s::Array # Scattering opacity
    energydep # Deposited energy in a timestep
    emittedenergy # Emitted energy during sourcing in a timestep
    totalenergy # Total energy emitted up to current time
    totalenergydep # Total energy deposited up to current time
    lostenergy # Energy of particles that left the problem domain 
    bee::Array # Heat Capacity
    radsource::Array # Radiation source term
    radenergydens::Array # Radiation energy density
    matenergydens::Array # Material energy density
    temp_saved # Temperature saved values 
    energyincrease_saved # Energy increase saved values
    radenergy_saved # Radiation energy density saved values
    matenergy_saved # Material energy density saved values
    energyscales # Scale factors for energy
    distancescale # Scale factor for distance
    centers # Center points on mesh
    nodes # Nodal points on mesh
    dx # Cell widths in x-direction
    dy # Cell widths in y-direction
    perimeter # Perimeter of the mesh (Tuple for 1D, Array for 2D)
    Ncells # Number of cells (Int or Tuple)
    end

    function mesh_generation(inputs)
        """
        Function to generate mesh quantities
        Parameters: 
        inputs: Dict - Dictionary of input parameters
        Returns:
        mesh: MeshStruct - Mesh object
        """

        print("Generating mesh quantities... \n")

        geometry = uppercase(inputs["GEOMETRY"])
        meshtype = uppercase(inputs["MESHTYPE"])
        precision = inputs["PRECISION"]

        if geometry == "1D"
            if meshtype == "UNIFORM"
                xsize = parse(precision, inputs["XSIZE"])
                dx_temp = parse(precision, inputs["DX"])
                Ncells =  Utilities.tointeger(xsize/dx_temp)
                dx = fill(dx_temp, Ncells)
                dy = 1
                nodes = LinRange(0, xsize, Ncells+1)
                centers =  LinRange(dx_temp/2, xsize - dx_temp/2, Ncells)
                perimeter = precision.((0,0))
            elseif meshtype == "NONUNIFORM"
                nodes = inputs["MESHNODES"]
                centers = precision.((nodes[1:end-1] + nodes[2:end]) / 2)
                dx = precision.(nodes[2:end] - nodes[1:end-1])
                dy = precision(1)
                nodes = precision.(nodes)
                Ncells = Utilities.tointeger(precision(length(dx)))
                perimeter = precision.((0,0))
            end
            temp_saved = Vector{Vector{precision}}()
            energyincrease_saved = Vector{Vector{precision}}()
            radenergy_saved = Vector{Vector{precision}}()
            matenergy_saved = Vector{Vector{precision}}()
        elseif geometry == "2D"
            if meshtype == "UNIFORM"
                xsize = parse(precision, inputs["XSIZE"])
                ysize = parse(precision, inputs["YSIZE"])
                dx_temp = parse(precision, inputs["DX"])
                dy_temp = parse(precision, inputs["DY"])
                Ncells = (Utilities.tointeger(xsize/dx_temp), Utilities.tointeger(ysize/dx_temp)) # Nx and Ny
                dx = fill(dx_temp, Ncells[1])
                dy = fill(dy_temp, Ncells[2])
                xnodes = LinRange(0, xsize, Ncells[1]+1)
                xcenters =  LinRange(dx_temp/2, xsize - dx_temp/2, Ncells[1])
                ynodes = LinRange(0, ysize, Ncells[2]+1)
                ycenters =  LinRange(dy_temp/2, ysize - dy_temp/2, Ncells[2])
                centers = (xcenters, ycenters)
                nodes = (xnodes, ynodes)
                perimeter = (zeros(precision, Ncells[1]), zeros(precision, Ncells[1]), zeros(precision, Ncells[2]), zeros(precision, Ncells[2])) # Bottom, Top, Left, Right
            elseif meshtype == "NONUNIFORM"
                xnodes = inputs["XMESHNODES"]
                ynodes = inputs["YMESHNODES"]
                xcenters = precision.((xnodes[1:end-1] + xnodes[2:end]) / 2)
                ycenters = precision.((ynodes[1:end-1] + ynodes[2:end]) / 2)
                centers = (xcenters, ycenters)
                dx = precision.(xnodes[2:end] - xnodes[1:end-1])
                dy = precision.(ynodes[2:end] - ynodes[1:end-1])
                xnodes = precision.(xnodes)
                ynodes = precision.(ynodes)
                nodes = (xnodes, ynodes)
                Ncells = (length(dx), length(dy))
                perimeter = (zeros(precision, Ncells[1]), zeros(precision, Ncells[1]), zeros(precision, Ncells[2]), zeros(precision, Ncells[2]))
            end
            temp_saved = Vector{Matrix{precision}}()
            energyincrease_saved = Vector{Matrix{precision}}()
            radenergy_saved = Vector{Matrix{precision}}()
            matenergy_saved = Vector{Matrix{precision}}()
        end

        T_init = parse(precision, inputs["T_INIT"])
        T_surface = surface_definer(geometry, inputs["T_SURFACE_REGS"], inputs["T_SURFACE_VALS"], perimeter, nodes, Ncells, precision) # Surface Temperature (keV)

        temp = fill(precision(T_init), Ncells) # Temperature (keV)
        push!(temp_saved, copy(temp)) # Save initial temperature distribution
        fleck = zeros(precision, Ncells) # Fleck factor
        beta = ones(precision, Ncells) # Beta factor

        sigma_a = region_joiner(geometry, inputs["SIGMA_A_REGS"], inputs["SIGMA_A_VALS"], nodes, Ncells, precision) # Absorption opacity
        sigma_a_powers = region_joiner(geometry, inputs["SIGMA_A_REGS"], inputs["SIGMA_A_POWERS"], nodes, Ncells, precision) # Absorption Opacity Temperature Power Laws
    
        sigma_s = region_joiner(geometry, inputs["SIGMA_S_REGS"], inputs["SIGMA_S_VALS"], nodes, Ncells, precision) # Scattering opacity
        sigma_s_powers = region_joiner(geometry, inputs["SIGMA_S_REGS"], inputs["SIGMA_S_POWERS"], nodes, Ncells, precision) # Scattering Opacity Temperature Power Laws
       
        energyscales = sort(inputs["ENERGYSCALES"], rev=true) # Scale factors for energy, sorted in descending order
        distancescale = parse(precision, inputs["DISTANCESCALE"]) # Scale factor for distance

        sigma_a = sigma_a / distancescale # Scale absorption opacity by distance scale factor
        sigma_s = sigma_s / distancescale # Scale scattering opacity by distance scale factor

        if geometry == "1D"
            sigma_a = hcat(sigma_a, sigma_a, sigma_a_powers) # Create a stacked array to store the full functional absorption opacity, its constants, and its powers
            sigma_s = hcat(sigma_s, sigma_s, sigma_s_powers) # Create a stacked array to store the full functional scattering opacity, its constants, and its powers
        elseif geometry == "2D"
            sigma_a = cat(sigma_a, sigma_a, sigma_a_powers, dims=3) # Create a stacked array to store the full functional absorption opacity, its constants, and its powers
            sigma_s = cat(sigma_s, sigma_s, sigma_s_powers, dims=3) # Create a stacked array to store the full functional scattering opacity, its constants, and its powers
        end

        sigma = sigma_a + sigma_s # Total opacity
        radsource = region_joiner(geometry, inputs["RADSOURCE_REGS"], inputs["RADSOURCE_VALS"], nodes, Ncells, precision) # Radiation source term

        if length(energyscales) == 1
            energyscales = parse(precision, inputs["ENERGYSCALES"][1])
        end

        if geometry == "1D"
            energydep = zeros(precision, (Ncells, length(energyscales))) # Energy depositions in a timestep, separate arrays for each energy scale 
            emittedenergy = zeros(precision, (Ncells, length(energyscales))) # Emitted energy during sourcing in a timestep
        elseif geometry == "2D"
            energydep = zeros(precision, (Ncells[1], Ncells[2], length(energyscales))) # Energy depositions in a timestep, separate arrays for each energy scale
            emittedenergy = zeros(precision, (Ncells[1], Ncells[2], length(energyscales))) # Emitted energy during sourcing in a timestep
        end
        
        totalenergy = precision(0.0)
        totalenergydep = precision(0.0)
        lostenergy = precision(0.0)

        bee = region_joiner(geometry, inputs["BEE_REGS"], inputs["BEE_VALS"], nodes, Ncells, precision) # Heat Capacity

        radenergydens = zeros(precision, Ncells) # Radiation energy density
        matenergydens = zeros(precision, Ncells) # Material energy density

        mesh = MeshStruct(temp, T_surface, fleck, beta, sigma, sigma_a, sigma_s, energydep, emittedenergy, totalenergy, totalenergydep, lostenergy, bee, radsource, radenergydens, 
                matenergydens, temp_saved, energyincrease_saved, radenergy_saved, matenergy_saved, energyscales, distancescale, centers, nodes, dx, dy, perimeter, Ncells)

        return mesh
    end

    function region_joiner(geometry, regions, values, nodes, Ncells, precision)
        """ Function to join regions and values for mesh generation
            Parameters:
            geometry: String - Geometry of the problem
            regions: Array - Array of regions
            values: Array - Array of values
            nodes: Array - Array of nodal points
            Ncells: Int - Number of cells
            precision: Type - Precision of the values
            Returns:
            joined_array: Array - Array of joined values
        """

        joined_array = zeros(precision, Ncells)

        region_index = 1
        if geometry == "1D"
            if length(regions) <= 1
                joined_array = fill(parse(precision, values[1]), Ncells)
            else
                for i in eachindex(joined_array)
                    if regions[region_index] >= nodes[i+1]
                        joined_array[i] = values[region_index]
                    else
                        region_index += 1
                        joined_array[i] = values[region_index]
                    end
                end
            end
        elseif geometry == "2D"
            if length(regions[1]) <= 1
                joined_array = fill(parse(precision, values[1]), Ncells)
            else
                for region_index in eachindex(regions[1])
                        xstart, xend = regions[1][region_index][1]
                        ystart, yend = regions[1][region_index][2]
                        xstartindex = findlast(x -> x <= xstart, nodes[1])
                        xendindex = findlast(x -> x <= xend, nodes[1])-1
                        ystartindex = findlast(y -> y <= ystart, nodes[2])
                        yendindex = findlast(y -> y <= yend, nodes[2])-1
                        joined_array[xstartindex:xendindex, ystartindex:yendindex] .= values[region_index]
                end
            end
        end
        return joined_array
    end
       
    function surface_definer(geometry, regions, values, perimeter, nodes, Ncells, precision)
        """ Function to define the surface temperature
            Parameters:
            geometry: String - Geometry of the problem
            regions: Array - Array of regions
            values: Array - Array of values
            perimeter: Array - Array of perimeter values
            nodes: Array - Array of nodal points
            Ncells: Int - Number of cells
            precision: Type - Precision of the values
            Returns:
            surface: Array - Array of surface values
        """
        if geometry == "1D"
            surface = (values[1], values[2])
        elseif geometry == "2D"
            surface = [[], [], [], []]
            surfacevals = [[], [], [], []]
            # Format values for region joiner function
            for i in 1:4
                if length(values[1][i]) <= 1
                    surfacevals[i] = [string(values[1][i])]
                else
                    surfacevals[i] = collect(values[1][i])
                end
            end
            # Bottom Boundary
            surface[1] = region_joiner("1D", collect(regions[1][1]), surfacevals[1], nodes[1], Ncells[1], precision)
            # Top Boundary
            surface[2] = region_joiner("1D", collect(regions[1][2]), surfacevals[2], nodes[1], Ncells[1], precision) 
            # Left Boundary
            surface[3] = region_joiner("1D", collect(regions[1][3]), surfacevals[3], nodes[2], Ncells[2], precision)
            # Right Boundary
            surface[4] = region_joiner("1D", collect(regions[1][4]), surfacevals[4], nodes[2], Ncells[2], precision)
        end

        return surface 
    end

end






