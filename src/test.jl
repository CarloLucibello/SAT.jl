include("SAT.jl")
using .SAT

N = 100
planted = [rand(1:N, N).*rand([-1,1], N) for _ in 1:10]
println(planted)
CNF = randomcnf(N=N, α=1.2)
println(CNF.clauses)


sol = solve(CNF, method=:converge)
println(sol)





