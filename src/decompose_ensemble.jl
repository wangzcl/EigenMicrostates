"""
    EigenEnsemble(A)

Decompose eigen ensembles (eigen microstates and their evolution) from an statistical ensemble array `A`. 

The last dimension of `A` should be the time dimension.
The statistical ensemble array `A` does not need to be normalized.

See https://doi.org/10.1088/1572-9494/abf127 for more details.

If `S::EigenEnsemble` represents eigen ensmebles of a `N`-dimensional system 
(e.g. a canonical ensmeble of an `N`-dimensional Ising model), then:

`S.U` is a `N+1`-dimensional array, `S.U[:, ..., :, i]` is the `i`-th eigen microstate;

`S.sigma` is a vector, representing propability amplitudes of eigen microstates, `sum(S.sigma.^2) == 1.0`;

`S.V` is a matrix, each column representing the dynamic evolution of the corresponding eigen microstate.

Eigen ensembles are arranged in descending order of propability amplitudes.

Use `ndims(S)` to get the dimension of the system under study (`N`).
"""
struct EigenEnsemble{T,Tr<:AbstractFloat,TU<:AbstractArray{T},Ts<:AbstractVector{Tr},TV<:Union{AbstractArray{T},Adjoint{T}},N}
    U::TU
    sigma::Ts
    V::TV
    EigenEnsemble(U::AbstractArray{T}, sigma::AbstractVector{Tr}, V::Union{AbstractArray{T},Adjoint{T}}) where {T,Tr} =
        new{T,Tr,typeof(U),typeof(sigma),typeof(V),ndims(U)-1}(U, sigma, V)
end


Base.ndims(::EigenEnsemble{T,Tr,TU,Ts,TV,N}) where {T,Tr,TU,Ts,TV,N} = N


function EigenEnsemble(A::AbstractArray)
    system_size = size(A)[begin:end-1]
    A = reshape_as_matrix(A)
    U, sigma, V = svd(A)
    U = reshape(U, system_size..., :)
    normalize!(sigma)
    return EigenEnsemble(U, sigma, V)
end


"""
    amplitudes_ensemble(A) -> sigma

Probability amplitude vector `sigma` of eigen ensembles of ensemble array `A`.

`A` does not need to be normalized.
"""
amplitudes_ensemble(A::AbstractArray) = normalize!(svdvals(reshape_as_matrix(A)))


"""
    weights_ensemble(A) -> weights

Probability ("weights") of eigen ensembles of ensmeble array `A`.
`weights = sigma.^2`.

`A` does not need to be normalized.
"""
weights_ensemble(A::AbstractArray) = amplitudes_ensemble(A) .^ 2
