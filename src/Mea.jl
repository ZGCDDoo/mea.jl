module Mea

include("Green.jl")
include("Periodize.jl")
include("Sigmadc.jl")

#using Mea.Green, Mea.Sigmadc, Mea.Periodize
export

# green.jl

# periodize.jl
ModelVector, buildmodelvec, calcdos, calcdos2,

#sigmadc.jl
calc_sigmadc


end # module
