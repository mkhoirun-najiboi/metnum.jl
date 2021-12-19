"""
    mulaipraktikum()
untuk membuka folder praktikum via IJulia.
"""
function mulaipraktikum()
    repo_directory = joinpath(@__DIR__,"..")
    path = joinpath(repo_directory, "notebookpraktikum")
    IJulia.notebook(;dir=path)
end
