using SAT
using Base.Test

cnf = CNF([[1, -5, 4], [-1, 5, 3, 4], [-3, -4]])
@test cnf.M == 3
@test cnf.N == 5
σ = solve(cnf);
@test energy(cnf, σ) == 0

cnf = randomcnf(N=2000, α=9.4, k=4, seed=19)
σ = solve(cnf, rstep=0.001, maxiters=1000);
@test energy(cnf, σ) == 0


# E, σ = solve(N=10000, α=9.6, k=4, seed_cnf=19, rstep=0.0002, maxiters=1000);
# @test E == 0

# DECIMATION NOT WORKING
# E, σ = KSAT.solve(method=:decimation, N=10000,α=9.6, k=4, seed_cnf=19,
#         r=0.02, maxiters=1000);
# @test E == 0
