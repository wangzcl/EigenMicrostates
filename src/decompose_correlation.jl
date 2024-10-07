"""
    EigenMicrostate(C, system_size)

Eigen microstates and weights of correlation matrix `C`.

`system_size` is the size of the system under study, used to reshape the eigen microstates.
`C` is transfered to a symmetric matrix first, in which only upper triangular part of `C` is used.

See https://doi.org/10.1007/s11433-018-9353-x for more details.

If `S::EigenMicrostate` represents eigen microstates of a `N`-dimensional system
(e.g. a canonical ensmeble of an `N`-dimensional Ising model), then:

`S.U` is a `N`-dimensional array, `S.U[:, ..., :, i]` is the `i`-th eigen microstate;

`S.weights` is a vector, representing propability of eigen microstates, `sum(S.weights) == 1.0`.

Eigen microstates are arranged in descending order of propability.

Use `ndims(S)` to get the dimension of the system under study (`N`).
"""
struct EigenMicrostate{T,Tr<:AbstractFloat,TU<:AbstractArray{T},Tw<:AbstractVector{Tr},N}
    U::TU
    weights::Tw
    EigenMicrostate(U::AbstractArray{T}, weights::AbstractVector{Tr}) where {T,Tr} =
        new{T,Tr,typeof(U),typeof(weights),ndims(U) - 1}(U, weights)
end

Base.ndims(::EigenMicrostate{T,Tr,TU,Tw,N}) where {T,Tr,TU,Tw,N} = N

function EigenMicrostate(C::AbstractMatrix, system_size::Dims{N}) where {N}
    C = Symmetric(C)
    e = eigen(C)
    weights = normalize!(reverse!(e.values), 1)
    U = reverse!(e.vectors; dims=2)
    U = reshape(U, system_size..., :)
    return EigenMicrostate(U, weights)
end


"""
    weights_correlation(C::AbstractMatrix) -> weights
Probability ("weights") of eigen microstates via decomposing correlation matrix `C`.
"""
weights_correlation(C::AbstractMatrix) = normalize!(reverse!(eigvals(Symmetric(C))), 1)


"""
    amplitudes_correlation(C::AbstractMatrix) -> amplitudes
Probability amplitudes of eigen microstates via decomposing correlation matrix `C`.
"""
function amplitudes_correlation(C::AbstractMatrix)
    w = weights_correlation(C)
    w[w.<0] .= zero(eltype(w))
    amp = sqrt.(w)
    return amp
end

"""
    order_param_correlation(C::AbstractMatrix) -> order_param

Order parameter as maximum eigenvalue of correlation matrix `C`.
"""
order_param_correlation(C::AbstractMatrix) = eigmax(Symmetric(C)) / tr(C)
