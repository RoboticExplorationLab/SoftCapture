

function uniquify(filename)
    while isfile(filename)
        split_filename = split(filename, ".")
        if occursin("(", split_filename[1])
            # update the value inside the parentheses by 1
            split_filename[1] = split(split_filename[1], " ")[1]*" ($(parse(Int, split(split_filename[1], " ")[2][2:end-1])+1)"*")"
        else
            # add a new parentheses with value 1
            split_filename[1] = split_filename[1]*" (1)"
        end
        # recombine the filename
        filename = join(split_filename, ".")
    end
    # # recursively call uniquify until unique filename is found
    # uniquify(filename)
    return filename
end

function mat_from_vec(X::Vector{Vector{Float64}})::Matrix
    # convert a vector of vectors to a matrix
    Xm = hcat(X...)
    return Xm
end

function vec_from_mat(Xm::Matrix)::Vector{Vector{Float64}}
    # convert a matrix into a vector of vectors 
    X = [Xm[:,i] for i = 1:size(Xm,2)]
    return X 
end


function has_nonzero_elements(matrix::AbstractMatrix, axis::Symbol)
    if axis == :row
        return all(x -> any(≠(0), x), eachrow(matrix))
    elseif axis == :column
        return all(x -> any(≠(0), x), eachcol(matrix))
    else
        throw(ArgumentError("Invalid axis argument. Use :row or :column."))
    end
end
