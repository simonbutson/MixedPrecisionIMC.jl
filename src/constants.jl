# Code to store various constant values
# Simon Butson

module Constants

function set_constants(inputs)
    precision = inputs["PRECISION"]
    # Physical Constants
    global phys_c = parse(precision,inputs["PHYS_C"]) #299.70 # speed of light [cm/shake]
    global phys_a = parse(precision,inputs["PHYS_A"]) #0.01372016 # radiation constant [Jrk/cm^3/K^4]

    global alpha = parse(precision,inputs["ALPHA"])
end

export alpha, phys_c, phys_a

end