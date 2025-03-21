# IMC Update Function
# Simon Butson

#print("phys_a is ", a, "\n")

module Update
using ..Constants

    function update(mesh, simvars)
        """
        Function to update the mesh parameters
        Parameters:
        mesh: Mesh - Mesh object
        simvars: SimVars - Simulation variables object
        Returns:
        None
        """
        precision = simvars.precision
        
        mesh.beta = (4 * phys_a * mesh.temp.^3) ./ mesh.bee

        if simvars.geometry == "1D"
            mesh.sigma_a[:,1] = mesh.sigma_a[:,2].*mesh.temp.^(mesh.sigma_a[:,3])
            #mesh.sigma_a[:,1] = mesh.sigma_a[:,2] ./ mesh.temp ./ mesh.temp ./ mesh.temp
            mesh.sigma_s[:,1] = mesh.sigma_s[:,2].*mesh.temp.^(mesh.sigma_s[:,3])
            sigma_a = mesh.sigma_a[:,1]
        elseif simvars.geometry == "2D"
            mesh.sigma_a[:,:,1] = mesh.sigma_a[:,:,2].*mesh.temp.^(mesh.sigma_a[:,:,3])
            mesh.sigma_s[:,:,1] = mesh.sigma_s[:,:,2].*mesh.temp.^(mesh.sigma_s[:,:,3])
            sigma_a = mesh.sigma_a[:,:,1]
        end

        #mesh.bee = precision.(4 * phys_a * mesh.temp.^3)

        #mesh.beta = ones(precision, mesh.Ncells)

        #print("sigma_a = ", sigma_a, "")

        mesh.fleck = precision.(1.0 ./ (1.0 .+ mesh.distancescale * alpha .* mesh.beta * phys_c * simvars.dt .* sigma_a))
        #print("Fleck factor = ", mesh.fleck, "\n")
    end
end