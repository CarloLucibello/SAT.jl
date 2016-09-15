# SAT
Heuristic algorithms based on message passing for boolean satisfaction problems (SAT).
## Installation
```
Pkg.clone("...")
```
## Usage
`solve` is the main function, and a solver `method` can be chosen
beetween `:reinforcement` (default),  `:decimation`

### Reinforcement
Solve random instance with BP inspired procedures.

`r` is the initial value of the reinforcement parameter (`r=0.` default).
`rstep` determines its moltiplicative increment.
```julia
E, σ = KSAT.solve(N=10000, α=9.6, k=4, seed_cnf=19, rstep=0.0002, maxiters=1000);
```

Read file in CNF format and solve with BP + reinforcement
```julia
E, σ = KSAT.solve("file.cnf", rstep=0.01, maxiters=1000);
```

If having errors, try to reduce `reinf_step`.

### Decimation
**NOT WORKING AT THE MOMENT**  
After each convergence of the BP algorithm the `r*N` most biased variables are fixed.
```julia
E, σ = KSAT.solve(method=:decimation, N=10000,α=9.6, k=4, seed_cnf=19, r=0.02, maxiters=1000);
