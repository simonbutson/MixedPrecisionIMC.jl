using MixedPrecisionIMC
using Documenter

DocMeta.setdocmeta!(MixedPrecisionIMC, :DocTestSetup, :(using MixedPrecisionIMC); recursive=true)

makedocs(;
    modules=[MixedPrecisionIMC],
    authors="simonbutson <butsons@oregonstate.edu> and contributors",
    sitename="MixedPrecisionIMC.jl",
    format=Documenter.HTML(;
        canonical="https://"simonbutson".github.io/MixedPrecisionIMC.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/"simonbutson"/MixedPrecisionIMC.jl",
    devbranch="main",
)
