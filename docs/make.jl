using JSMDUtils
using Documenter

makedocs(;
    modules=[JSMDUtils],
    authors="Andrea Pasquale <andrea.pasquale@polimi.it> and Michele Ceresoli <michele.ceresoli@polimi.it>",
    sitename="JSMDUtils.jl",
    pages=["Home" => "index.md", 
            "Utilities" => "utils.md"
        ],
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/JSMDUtils.jl", branch="gh-pages"
)