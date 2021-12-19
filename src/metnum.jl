"""
    metnum.jl

Modul ini berisi metode untuk menyelesaikan

1. Akar pertaksamaan tak-linear: `bisection`, `regulaFalsi`, `fixpoint`, `newtonRaphson`, `secant`
2. SPL langsung: `backsub`, `elimGaussNonPivoting`, `elimGaussWithPivoting`, `LUtanpaP`, `LUdenganP`
3. SPL iteratif: `jacobi`, `gaussSeidel`, `rekons` ,`conGrad`
4. Interpolas: `lagrange`, `newton`
5. Regresi: `reglin`, `regpower`, `regpoly`
6. Turunan: `bedaPusat`, `richardson`
7. Integral: `komptrap`, `kompsimp`, `rekursif`, `romberg`, `adaptif`
8. PDB dengan MNA: `euler`, `heun`, `taylor`, `rungekutta`, `rungekuttasistem`, `rkf45`
9. PDB dengan MNB: `linearshooting`, `findiff`
10. PDP: `gelombang`, `panas`, `dirichlet`
"""
module metnum

using LinearAlgebra
using Statistics
using Dates
using IJulia

export regpoly, regpower, regulaFalsi, rekons, rekursif, richardson, rkf45
export romberg, rungekutta, rungekuttasistem, secant
export taylor, adaptif, backsub, bedaPusat, bisection
export conGrad, dirichlet, elimGaussNonPivoting, elimGaussWithPivoting
export euler, findiff, fixpoint, gaussSeidel, gelombang, heun, jacobi
export kompsimp, komptrap, lagrange, linearshooting, LUdenganP, LUtanpaP
export newton, newtonRaphson, panas, reglin, mulaipraktikum

include("regpoly.jl")
include("regpower.jl")
include("regulaFalsi.jl")
include("rekons.jl")
include("rekursif.jl")
include("richardson.jl")
include("rkf45.jl")
include("romberg.jl")
include("rungekutta.jl")
include("rungekuttasistem.jl")
include("secant.jl")
include("taylor.jl")
include("adaptif.jl")
include("backsub.jl")
include("bedaPusat.jl")
include("bisection.jl")
include("conGrad.jl")
include("dirichlet.jl")
include("elimGaussNonPivoting.jl")
include("elimGaussWithPivoting.jl")
include("euler.jl")
include("findiff.jl")
include("fixpoint.jl")
include("gaussSeidel.jl")
include("gelombang.jl")
include("heun.jl")
include("jacobi.jl")
include("kompsimp.jl")
include("komptrap.jl")
include("lagrange.jl")
include("linearshooting.jl")
include("LUdenganP.jl")
include("LUtanpaP.jl")
include("newton.jl")
include("newtonRaphson.jl")
include("panas.jl")
include("reglin.jl")


function mulaipraktikum()
    using Dates
    repo_directory = joinpath(@__DIR__,".julia","packages","metnum")
    folders = readdir(repo_directory)
    unix_time = Array{Any}(undef,length(folders))
    k=1
    for folder in folders
        date_time = unix2datetime(mtime(joinpath(repo_directory, folder)))
        unix_time[k] = Int(floor(datetime2unix(date_time)))
        k+=1
    end
    path = joinpath(repo_directory, folders[argmax(unix_time)],"notebookpraktikum")
    IJulia.notebook(;dir=path)
end


end
