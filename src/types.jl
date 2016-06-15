type CSOut{T<:AbstractFloat}
    N::Int
    M::Int
    fmin::T
    x::Vector{T}
    z::Vector{T}    
    t::T
    status::Symbol
end

type RunOut{T<:AbstractFloat}
    res::CsensLP.CSOut{T}
    y::Vector{T}
    A::Matrix{T}
    x::Vector{T}
end

type ResOut{T<:AbstractFloat}
    vecene::Vector{T}
    vecdiff::Vector{T}
    vecspar::Vector{Int}
end
