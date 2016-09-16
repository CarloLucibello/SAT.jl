using Documenter
# include("../src/LightGraphs.jl")
using SAT

# index is equal to the README for the time being
cp(normpath(@__FILE__, "../../README.md"), normpath(@__FILE__, "../src/index.md"); remove_destination=true)

makedocs(modules=[SAT], doctest = false)

rm(normpath(@__FILE__, "../src/index.md"))

deploydocs(
    deps   = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    repo   = "github.com/CarloLucibello/SAT.jl.git",
    julia  = "release"
)
