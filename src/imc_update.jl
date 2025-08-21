# IMC Update Function
# Simon Butson

#print("phys_a is ", a, "\n")

module Update

include("imc_utilities.jl")
using .Utilities
using ..Constants

    function update(inputs, mesh, simvars)
        """
        Function to update the mesh parameters
        Parameters:
        mesh: Mesh - Mesh object
        simvars: SimVars - Simulation variables object
        Returns:
        None
        """
        precision = simvars.precision
        
        if inputs["LINEARIZED"] == "TRUE"
            mesh.bee .= 4 * phys_a * mesh.temp.^3
            mesh.beta = ones(precision, mesh.Ncells)
        else
            mesh.beta = (4 * phys_a * mesh.temp.^3) ./ mesh.bee
        end

        if simvars.geometry == "1D"            
            mesh.sigma_a[:,1] = mesh.sigma_a[:,2].*mesh.temp.^(mesh.sigma_a[:,3])
            if uppercase(inputs["NAME"]) == "MARSHAK WAVE"
                mesh.sigma_a[:,1] = mesh.sigma_a[:,2] ./ mesh.temp ./ mesh.temp ./ mesh.temp
            end
            mesh.sigma_s[:,1] = mesh.sigma_s[:,2].*mesh.temp.^(mesh.sigma_s[:,3])
            sigma_a = mesh.sigma_a[:,1]

            for i in 1:mesh.Ncells
                mesh.fleck[i] = precision(1.0 / (1.0 + Utilities.sorter([mesh.distancescale, alpha, mesh.beta[i], phys_c, simvars.dt, sigma_a[i]], [1] , precision)[1]))
                if mesh.fleck[i] == precision(1.0)
                    print("Fleck factor is invalid at cell ", i, " with values: ", mesh.distancescale, " ", alpha, " ", mesh.beta[i], " ", phys_c, " ", simvars.dt, " ", sigma_a[i], " Temperature: ", mesh.temp[i], "\n")
                end
            end

        elseif simvars.geometry == "2D"
            mesh.sigma_a[:,:,1] = mesh.sigma_a[:,:,2].*mesh.temp.^(mesh.sigma_a[:,:,3])
            mesh.sigma_s[:,:,1] = mesh.sigma_s[:,:,2].*mesh.temp.^(mesh.sigma_s[:,:,3])
            sigma_a = mesh.sigma_a[:,:,1]

            for i in 1:mesh.Ncells[1]
                for j in 1:mesh.Ncells[2]
                    mesh.fleck[i,j] = precision(1.0 / (1.0 + Utilities.sorter([mesh.distancescale, alpha, mesh.beta[i,j], phys_c, simvars.dt, sigma_a[i,j]], [1] , precision)[1]))
                end
            end
            #mesh.fleck = precision.(1.0 ./ (1.0 .+ mesh.distancescale * alpha .* mesh.beta * phys_c * simvars.dt .* sigma_a))
        end

        #mesh.bee = precision.(4 * phys_a * mesh.temp.^3)

       

        #print("sigma_a = ", sigma_a, "")

        

        #print(Utilities.sorter([mesh.distancescale, alpha, mesh.beta[1], phys_c, simvars.dt, sigma_a[1]], [1] , precision)[1], "\n")

        #mesh.fleck = precision.(1.0 ./ (1.0 .+ mesh.distancescale * alpha .* mesh.beta * phys_c * simvars.dt .* sigma_a))
        #print("Fleck factor = ", mesh.fleck, "\n")
    end
end