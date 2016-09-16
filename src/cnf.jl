"""
A type representing a conjunctive normal form.
"""
type CNF
    N::Int
    M::Int
    clauses::Vector{Vector{Int}}
end

function CNF(clauses::Vector{Vector{Int}})
    M = length(clauses)
    N = maximum(maximum.(abs.(clauses)))
    return CNF(N, M, clauses)
end

"""
    randomcnf(; N=100, k=3, α=0.1, seed=-1, planted = Vector{Vector{Int}}())

Generates a random instance of the k-SAT problem, with `N` variables
and `αN` clauses.

Any configuration in `planted` is guaranteed to be a solution of the problem.
"""
function randomcnf(; N::Int = 100, k::Int = 3, α::Float64 = 0.1, seed::Int=-1,
                    planted = Vector{Vector{Int}}())
    seed > 0 && srand(seed)
    M = round(Int, N*α)
    clauses = Vector{Vector{Int}}()
    for p in planted
        @assert length(p) == N   "Wrong size for planted configurations ($N != $(lenght(p)) )"
    end
    for a=1:M
        while true
            c = rand(1:N, k)
            length(union(c)) != k && continue
            c = c .* rand([-1,1], k)

            # reject if not satisfies the planted solutions
            sat = Bool[any(i -> i>0, sol[abs.(c)] .* c) for sol in planted]
            !all(sat) && continue

            push!(clauses, c)
            break
        end
    end
    return CNF(N, M, clauses)
end


"""
    readcnf(fname::String)

Reads a CNF from file `fname`.
"""
function readcnf(fname::String)
    f = open(fname, "r")
    head = split(readline(f))
    N, M = parse(Int, head[3]), parse(Int, head[4])
    clauses = Vector{Vector{Int}}()
    for i=1:M
        line = readline(f)
        c = [parse(Int64, e) for e in split(line)[1:end-1]]
        push!(clauses, c)
    end
    return CNF(N, M, clauses)
end


"""
    writecnf(fname::String, cnf::CNF)

Writes `cnf` to file `fname`.
"""
function writecnf(fname::String, cnf::CNF)
    f = open(fname, "w")
    println(f, "p cnf $(cnf.N) $(cnf.M)")
    for c in cnf.clauses
        for i in c
            print(f, i, " ")
        end
        print(f, "0\n")
    end
    close(f)
end
