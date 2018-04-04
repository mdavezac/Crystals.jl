using Documenter
using Crystals

makedocs(modules = [Crystals],
         clean = true,
         format = :html,
         sitename = "Crystals.jl",
         authors = "Mayeul d'Avezac",
         analytics = "UA-89508993-1",
         pages = Crystals._doc_pages())

deploydocs(
    repo = "github.com/mdavezac/Crystals.jl.git",
    target = "build",
    julia = "0.6",
    deps = nothing,
    make = nothing,
)
