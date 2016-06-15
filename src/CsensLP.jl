module CsensLP


using Gurobi
using JuMP
using Convex

# given y,A  find x s.t 1)  y = A * x , 2)  |x|_1 is minimal 
# gurobi and jump version should solve (it is not clear)
# given y,A,λ  find  argmin_x { ||y - A * x||_∞ , + λ |x|_1 }

export CSOut, RunOut, ResOut, cslconvex, cslgurobi, csljump, compute_roc, runanal, analdata, alphac

include("types.jl")
include("cs.jl")
include("runanal.jl")
include("l1curve.jl")

end #module end
