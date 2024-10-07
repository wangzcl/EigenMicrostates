module EigenMicrostates

using LinearAlgebra
include("utils.jl")
include("decompose_ensemble.jl")
include("decompose_correlation.jl")

export EigenEnsemble
export amplitudes_ensemble, weights_ensemble

export EigenMicrostate
export amplitudes_correlation, weights_correlation, order_param_correlation

end
