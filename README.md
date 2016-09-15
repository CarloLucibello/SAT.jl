# SAT

[![Build Status](https://travis-ci.org/CarloLucibello/SAT.jl.svg?branch=master)](https://travis-ci.org/CarloLucibello/SAT.jl)
[![codecov.io](http://codecov.io/github/CarloLucibello/SAT.jl/coverage.svg?branch=master)](http://codecov.io/github/CarloLucibello/SAT.jl?branch=master)

Heuristic algorithms based on message passing for solving large instances of boolean satisfaction problems (SAT).

## Installation
```
Pkg.clone("https://github.com/CarloLucibello/SAT.jl")
```

## Basic usage
```julia
using SAT
cnf = randomcnf(N=1000, k=3, α=0.5) # generate a random k-SAT instance
sol = solve(cnf)
```
## CNF
Formulas in conjunctive normal form (![CNF](https://en.wikipedia.org/wiki/Conjunctive_normal_form)) can be either read/written to files
```julia
cnf = readcnf("formula.cnf")
writecnf("formula.cnf", cnf)
```
or generated randomly from the k-SAT ensemble
```julia
cnf = randomcnf(N=1000, k=4, α=0.5, seed=17)
```

## Belief Propagation
Solve random instance with Belief Propagation (BP) inspired procedures.

### reinforcement
`r` is the initial value of the reinforcement parameter (`r=0.` default).
`rstep` determines its moltiplicative increment.
```julia
E, σ = solve(cnf, rstep=0.001, maxiters=1000);
```
If having errors or unable to find a solution, try to reduce `rstep`.
### decimation
TODO
