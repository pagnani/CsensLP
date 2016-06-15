function runanal(M,  N,  K, nsample)    
    res=CSOut[]
    myprocs = procs()
    nproc = length(myprocs)
    if nproc == 1
        for s=1:nsample
            y,A,x = creamatrix(M,N,K)
            _res = cslconvex(y,A,verbose=0)
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
    vecdiff = Float64[norm(res[i].x - res[i].res.x) for i=1:length(res)]
    vecspar = Int[sum(abs(res[i].res.x) .<= 1e-10) for i=1:length(res)] # 1e-10 is should be slightly larger the numerical precision of the algorithm ( O(1e-12) )
    ResOut(vecene, vecdiff,vecspar)
end


nothing
