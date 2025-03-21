# IMC Cleaning Routine
# Simon Butson

module Clean

function clean(particles)
    """ Remove particles that have exited the problem domain or have been absorbed 
    Parameters:
    particles: Array{Array{Float64,1},1} - Array of particle data
    Returns:
    None
    """
    for prtclindex in length(particles):-1:1 # Loop backwards
        #print(particles[prtclindex][:], "\n")
        if particles[prtclindex][8] == -1.0
            deleteat!(particles, prtclindex)
        end
    end
end


end