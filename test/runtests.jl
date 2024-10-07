using EigenMicrostates
using Test

@testset "EigenMicrostates.jl" begin
    ising_ensemble = [
        1 1 1 -1
        -1 1 1 -1
        1 -1 -1 -1
        1 -1 -1 1
        -1 -1 1 1
        -1 -1 1 -1
    ]
    N, M = size(ising_ensemble) # N = 6, M = 4
    S = EigenEnsemble(ising_ensemble)
    @test ndims(S) == 1
    @test size(S.U) == (N, M)
    @test size(S.sigma) == (M,)
    @test size(S.V) == (M, M)

    @test amplitudes_ensemble(ising_ensemble) ≈ S.sigma
    @test weights_ensemble(ising_ensemble) ≈ S.sigma .^ 2

    ising_correlation = ising_ensemble * ising_ensemble'
    C = EigenMicrostate(ising_correlation, (size(ising_ensemble, 1),))
    @test ndims(C) == 1
    @test size(C.U) == (N, N)
    @test size(C.weights[begin:begin-1+M]) == size(S.sigma .^ 2)

    @test weights_correlation(ising_correlation) ≈ C.weights
    @test amplitudes_correlation(ising_correlation).^2 ≈ C.weights

end
