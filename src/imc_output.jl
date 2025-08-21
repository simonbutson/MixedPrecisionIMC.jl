# IMC Output Files 
# Simon Butson

module Output
using Plots
using JLD2

function plotting(inputs, mesh, simvars)
    """ Function to plot the output of the IMC simulation
        Parameters:
        inputs: Dict - Dictionary of input parameters
        mesh: Mesh - Mesh object
        simvars: SimVars - Simulation variables object
        Returns:
        None
    """

    if inputs["PLOTVARS"] == "TEMPERATURE"
        p = plot(mesh.centers, mesh.temp, markershape=:circle, title=inputs["NAME"]*" "*string(inputs["PRECISION"])*" at t = "*string(round(simvars.t, sigdigits=3)), label="IMC")
        xlabel!("x [cm]")
        ylabel!("Temperature [keV]")
        if inputs["BENCHMARK"] == "TRUE"
            plot!(inputs["XBENCH"], inputs["YBENCH"], markershape=:cross, label="Benchmark")
        end
        display(p)
        
    elseif inputs["PLOTVARS"] == "RADENERGY"
        p = plot(mesh.centers, mesh.radenergydens, markershape=:circle, title=inputs["NAME"]*" "*string(inputs["PRECISION"])*" at t = "*string(round(simvars.t, sigdigits=3)), label="IMC")
        xlabel!("x [cm]")
        ylabel!("Radiation Energy Density [keV/cm^3]")
        if inputs["BENCHMARK"] == "TRUE"
            plot!(inputs["XBENCH"], inputs["YBENCH"], markershape=:cross, label="Benchmark")
        end
        display(p)
    elseif inputs["PLOTVARS"] == "MATENERGY"
        p = plot(mesh.centers, mesh.matenergydens, markershape=:circle, title=inputs["NAME"]*" "*string(inputs["PRECISION"])*" at t = "*string(round(simvars.t, sigdigits=3)), label="IMC")
        xlabel!("x [cm]")
        ylabel!("Material Energy Density [keV/cm^3]")
        if inputs["BENCHMARK"] == "TRUE"
            plot!(inputs["XBENCH"], inputs["YBENCH"], markershape=:cross, label="Benchmark")
        end
        display(p)
    elseif inputs["PLOTVARS"] == "HEATMAP"
        p = heatmap(mesh.centers[1], mesh.centers[2], mesh.temp', title=inputs["NAME"]*" "*string(inputs["PRECISION"])*" at t = "*string(round(simvars.t,sigdigits=3)),  colorbar_title=" \nTemperature [keV]", clims=(inputs["CLIMS"][1],inputs["CLIMS"][2]), right_margin=5Plots.mm)
        xlabel!("x [cm]")
        ylabel!("y [cm]")
        display(p)
    end

    if simvars.t == simvars.t_end 
        if inputs["SAVEVARS"] == "TRUE"

                        
            # Save the mesh and simulation variables to a JLD2 file
            filename = "..\\outputs\\" * inputs["NAME"] * "_" * string(inputs["PRECISION"]) * ".jld2"
            jldsave(filename; simvars.timesteps, mesh.centers, mesh.temp_saved, mesh.radenergy_saved, mesh.matenergy_saved)

    
            # open("..\\outputs\\" * inputs["NAME"] * "_" * string(inputs["PRECISION"]) * ".txt", "w") do io
            #     write(io, "Time: ", simvars.t, "\n")
            #     write(io, "Temperature: ", mesh.temp, "\n")
            #     write(io, "Radiation Energy Density: ", mesh.radenergydens, "\n")
            #     write(io, "Material Energy Density: ", mesh.matenergydens, "\n")
            # end
        end
        if inputs["SAVEFIG"] == "TRUE"
            savefig("..\\outputs\\" * inputs["NAME"] * "_" * inputs["PLOTVARS"] * "_" * string(inputs["PRECISION"]) * ".png")
        end
        if inputs["SAVEANIMATION"] == "TRUE"
            print(length(simvars.timesteps), " frames to generate \n")
            anim = @animate for i = 1:length(simvars.timesteps)
                print("Generating frame ", i, "\n")
                heatmap(mesh.centers[1], mesh.centers[2], mesh.temp_saved[i]', fill= true, title=inputs["NAME"]*" "*string(inputs["PRECISION"])*" at t = "*string(round(simvars.timesteps[i],sigdigits=3)), colorbar_title=" \nTemperature [keV]", clims=(0,0.5), right_margin=5Plots.mm)
                xlabel!("x [cm]")
                ylabel!("y [cm]")
            end
            gif(anim, "..\\outputs\\" * inputs["NAME"] * "_" * inputs["PLOTVARS"] * "_"*string(inputs["PRECISION"]) * ".gif", fps = 5)
        end
    end
end

end