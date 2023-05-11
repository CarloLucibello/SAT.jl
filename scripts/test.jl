using StochasticLocalization

N = 100
# planted = [rand(1:N, N).*rand([-1,1], N) for _ in 1:10]
# println(planted)
CNF = randomcnf(; N, α=3.5, k=3, seed=1)
println(CNF.clauses)

sol = solve(CNF, method=:reinforcement, r=0., rstep=0.01, 
            seed=-1, infotime=1, 
            altsolv=true,
            ϵ=1e-8,
            maxiters=100000);
# println(sol)





