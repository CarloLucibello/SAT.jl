module SAT
using Random
using ExtractMacro
using Printf

export solve, energy,
        CNF, randomcnf, readcnf, writecnf

include("cnf.jl")
include("belief_propagation.jl")
include("gradient_descent.jl")

end # module
