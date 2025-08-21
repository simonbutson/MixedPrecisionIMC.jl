# IMC Utility Functions
# Simon Butson

module Utilities

function tointeger(x::AbstractFloat)
    """ Function to convert a float to an integer of the same precision
        Parameters:
        x: AbstractFloat - Float to be converted
        Returns:
        Int - Integer value of x with the same precision
    """

    if typeof(x) == Float16
        return UInt16(x)
    elseif typeof(x) == Float32
        return Int32(x)
    else
        return Int64(x)
    end
end

function sorter(vals, scales, precision)
    """ Function to sort the values and change order of operations 
        to return the product of all the values
        Parameters:
        vals: Vector{Float} - Vector of values to be sorted
        scales: Vector{Float} - Vector of scale factors
        Returns:
        Tuple{Float, Float, Int} - Tuple containing the product of the values, scale factor, and index of the scale factor
    """

    #scales = sort(scales, rev=true) # Sort the scales in descending order

    for j in eachindex(scales)
        product = precision(1.0)
        sorted = sort(vcat(vals, scales[j]))

        for i in 1:div(length(sorted),2)
            product *= precision(sorted[i] * sorted[end-i+1])
        end
        if isodd(length(sorted))
            product *= precision(sorted[div(length(sorted),2)+1])
        end
        if !isinf(product) && !isnan(product)
            return (product, scales[j], j)
        end
    end

    print("Product not representable \n")
    print(vals, scales, "\n")
    sleep(1.0) # Pause to allow user to see the error
    return (precision(0.0), precision(0.0), 0) # Return a default value if no product is found
end


end