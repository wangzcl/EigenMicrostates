"""
    reshape_as_matrix(A)
Reshape an array `A` of dimension N (N>1) into a matrix, by collapsing
all the dimensions except the last one into a single dimension.
"""
function reshape_as_matrix(A::AbstractArray)
    @assert ndims(A) > 1
    return reshape(A, :, size(A)[end])
end