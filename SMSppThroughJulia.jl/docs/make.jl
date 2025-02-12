push!(LOAD_PATH,joinpath(@__DIR__, ".."))
using Documenter, SMSppThroughJulia

makedocs(
    modules = [SMSppThroughJulia],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Wim van Ackooij",
    sitename = "SMSppThroughJulia.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com//SMSppThroughJulia.jl.git",
)
