function cslconvex(y::Vector{Float64},A::Matrix{Float64};
                   verbose::Int=1, 
                   solver::Symbol=:gurobi #the other option are :glpk
                   )

    if solver == :gurobi
        strsolv = GurobiSolver(OutputFlag=verbose)
    elseif solver == :glpk
        strsolv = GLPKSolverMIP()
    else
        error("solver $solver seems to not being installed\n")
    end

    M,N = size(A)
    mytime = @elapsed begin
        xx = Convex.Variable(N)
        problem = minimize(norm(xx,1), A*xx == y)
        solve!(problem,strsolv)
    end
    println("problem solved in $mytime [s]")
    CSOut(M,N,problem.optval, vec(xx.value), zeros(N), 0.0, problem.status) 
end




function cslgurobi(y::Vector{Float64},A::Matrix{Float64};
                   ρ::Real=0.5 # uniform sparsity
                   )

    M,N = size(A)
    length(y) == M || error("number of rows of A = $M different from length(y) = ",length(y))

    if ρ > 0
        coeffopt = Float64[ρ * Float64(i>N) for i=1:2N+1]
    elseif ρ == 0
        coeffopt = zeros(Float64, 2N+1)
    else
        error("ρ = $ρ should be positive")
    end       

    coeffopt[end]=1.0
        
    env = Gurobi.Env()
    
    m = gurobi_model(env; sense=:minimize, f=coeffopt)
    
    for i=1:N
        add_constr!(m,[i+N,i],[1.0,-1.0],'>', 0.0) 
        add_constr!(m,[i+N,i],[1.0, 1.0],'>', 0.0)
    end
    myind = collect(1:N+1)
    myind[end] = 2N+1
    mycoe = zeros(N+1)
    
    for i=1:M
        ai = vec(A[i,:])
        for j=1:N
            mycoe[j]=-ai[j] 
        end
        mycoe[end] = 1.0
        add_constr!(m,myind, mycoe, '>', -y[i])
        mycoe[end] = -1.0
        add_constr!(m,myind, mycoe, '<', -y[i])
    end

    mytime = @elapsed optimize(m)
    println("elapsed time = $mytime")
    _x = get_solution(m)        
    CSOut(M,N, get_objval(m), _x[1:N],_x[N+1:2N],_x[2N+1], Gurobi.get_status(m))
    
end

function csljump(y::Vector{Float64},A::Matrix{Float64};
              ρ::Real=0.5 # uniform sparsity
              )
    
    M,N = size(A)
    length(y) == M || error("number of rows of A = $M different from length(y) = ",length(y))
    m = Model(solver = GurobiSolver())
    @variable(m,x[1:2N+1]) 
    for i=1:M
        ai = vec(A[i,:])
        @constraint(m, x[2N+1] >= (sum{ai[j]*x[j], j=1:N} - y[i]))
        @constraint(m,-x[2N+1] <= (sum{ai[j]*x[j], j=1:N} - y[i]))
    end
    for i=1:N
        @constraint(m,  x[i+N] >= x[i])
        @constraint(m, -x[i+N] <= x[i])
    end

    if ρ > 0.0
        @objective(m, Min, x[2N+1]+ρ*sum{x[i],i=N+1:2N})
    elseif ρ == 0
        @objective(m, Min, x[2N+1])
    else
        error("ρ = $ρ should be positive")
    end
    status= solve(m)
    
    res = getvalue(x)
    
    CSOut(M,N, getobjectivevalue(m), res[1:N],res[N+1:2N],res[2N+1], status)
    
end


function compute_roc(xtrue::Vector,x::Vector)
    idx = sortperm(abs(x),rev=true)
    tp = 0
    ctr = 0
    roc = zeros(length(idx))
    for i in idx
        ctr += 1
        if abs(xtrue[i]) > 0
            tp += 1
        end
        roc[ctr] = tp/ctr        
    end
    roc
end
