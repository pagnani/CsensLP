function runanal(M,  N,  K, nsample)    


    myprocs = procs()
    nproc = length(myprocs)
    if nproc == 1
        res=RunOut[]
        for s=1:nsample
            y,A,x = creamatrix(M,N,K)
            _res = RunOut(cslconvex(y,A,verbose=0), y, A,x)            
            push!(res,_res)
        end
    else
        res = @parallel vcat for s=1:nsample
            y,A,x = creamatrix(M,N,K)
            # y = rand(M)
            # x = rand(N)
            # A = rand(M,N)
            _res = cslconvex(y,A,verbose = 0)
            RunOut(_res, y, A,x)
        end
    end
    res
end

function creamatrix(M,N,K)
    A = randn(M,N)
    x = zeros(N)
    for i in randperm(N)[1:K]
        x[i] = randn()
    end
    y = A*x
    return y, A, x

end

function analdata(res::Vector)

    vecene  = Float64[norm(res[i].A * res[i].res.x - res[i].y) for i=1:length(res)]
    vecdiff = Float64[norm(res[i].x - res[i].res.x,1)/res[i].res.N for i=1:length(res)]
    vecspar = Int[sum(abs(res[i].res.x) .<= 1e-10) for i=1:length(res)] # 1e-10 is should be slightly larger the numerical precision of the algorithm ( O(1e-12) )
    ResOut(vecene, vecdiff,vecspar)
end

function runbatch(N::Int, rho0::Float64; nsample::Int=100,nstep::Int=7)
    
    alpha=linspace(rho0, alphac(rho0) + rho0, 7)
    vecM = round(Int, N * collect( alpha))
    K = round(Int, rho0 * N)

    res = Vector{Vector{RunOut}}()   
    out = Vector{ResOut}()
    for M in vecM
        _res = runanal(M, N, K, nsample) 
        _out = analdata(_res)
        push!(res,_res)
        push!(out,_out)
    end


    res,out
end


nothing
