module metnum

using LinearAlgebra
using Statistics

export regpoly, regpower, regulaFalsi, rekons, rekursif, richardson, rkf45
export romberg, rungekutta, rungekuttasistem, secant
export taylor, adaptif, backsub, bedaPusat, bisection
export conGrad, dirichlet, elimGaussNonPivoting, elimGaussWithPivoting
export euler, findiff, fixpoint, gaussSeidel, gelombang, heun, jacobi
export kompsimp, komptrap, lagrange, linearshooting, LUdenganP, LUtanpaP
export newton, newtonRaphson, panas, reglin

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


end
