using PyCall
@pyimport scipy.optimize as so
sqrt2 = sqrt(2.0)
H(x) = 0.5*erfc(x/sqrt2)
G(x) = sqrt(x/(2π)) * exp(-1/(2x)) - H(1/sqrt(x))


function l1curve(xinf,xsup; nstep=100)
    xx = linspace(xinf,xsup, nstep)
    for x in xx
        ρ =  2G(x)/(1+2G(x))
        α = 2(1-ρ)*H(1/sqrt(x)) + ρ
        println("G($x) = ",G(x)," ρ = ",ρ , " α = ", α)
    end
end

@vectorize_1arg Number alphac

function alphac(ρ)
    myrho(x) = 2G(x)/(1+2G(x)) - ρ
    
    if 0. > ρ > 1.
        error("ρ = $ρ should be in [0,1]")
    elseif ρ == 1.0
        x = 1.0
    elseif ρ == 0
        x = 0
    else
        ρ <= 0.5 ? (x0 = max(ρ,0.2)) : (x0 = min(ρ,0.9))
        x = so.newton(myrho, x0 , maxiter=1000000, tol=1e-15)
    end
    2(1-ρ)*H(1/sqrt(x)) + ρ

end



