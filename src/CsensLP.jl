module CsensLp

# given y,A  find x s.t 1)  y = A * x , 2)  |x|_1 is minimal 
# gurobi and jump version should solve (it is not clear)
# given y,A,λ  find  argmin_x { ||y - A * x||_∞ , + λ |x|_1 }

using Gurobi,JuMP,Convex

export CSOut, RunOut, ResOut, cslconvex, cslgurobi, csljump, compute_roc, runanal, analdata

include("types.jl")
include("cs.jl")
include("runanal.jl")




end #module end
