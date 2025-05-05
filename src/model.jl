nullfunc!(v, r)= (fill!(v, 0.0); return 0.0);

abstract type ODEModels end
mutable struct MosquitoModel <: ODEModels
    ndim::Int
    verbose::Bool
    islearned::MyTuple{Bool, Bool, Bool, Bool}
    islearning::MyTuple{Bool, Bool, Bool, Bool}
    funcs::MyTuple{Vector{Ψ1{F1}}, Vector{Ψ2{F2a,F2b,F2ab}}, Vector{Ψ2{F3a,F3b,F3ab}}, Vector{Ψ2{F4a,F4b,F4ab}}} where {F1,F2a,F2b,F2ab,F3a,F3b,F3ab,F4a,F4b,F4ab}
    coeffs::MyTuple{Vector{T}, Vector{T}, Vector{T}, Vector{T}} where T
    params::MyTuple{Vector{T}, Vector{T}, Vector{T}, Vector{T}} where T
    bfields::MyTuple{Function, Function, Function, Function}
    diffD::Union{T, Vector{T}} where T
    geom::Geometry
    velocity::Vector{T} where T
    position::Vector{T} where T
    bvec::Vector{T} where T
    uvec::Vector{T} where T
    MosquitoModel(;ndim=2, verbose=true,
        f1  =LaguerreFunction{Float64}(),
        f2a =LaguerreFunction{Float64}(),
        f2b =LaguerreFunction{Float64}(),
        f2ab=LegendreP{Float64}(),
        f3a =LaguerreFunction{Float64}(),
        f3b =LaguerreFunction{Float64}(),
        f3ab=LegendreP{Float64}(),
        f4a =LaguerreFunction{Float64}(),
        f4b =LaguerreFunction{Float64}(),
        f4ab=LegendreP{Float64}()
        ) = new( ndim, verbose, MyTuple(0,0,0,0), MyTuple(0,0,0,0),
                MyTuple(Vector{Ψ1{f1}}(undef,1),Vector{Ψ2{f2a,f2b,f2ab}}(undef,1),Vector{Ψ2{f3a,f3b,f3ab}}(undef,1),Vector{Ψ2{f4a,f4b,f4ab}}(undef,1)),
                MyTuple(zeros(1),zeros(1),zeros(1),zeros(1)),
                MyTuple(ones(1), ones(2), ones(2), ones(2)),
                MyTuple(nullfunc!,nullfunc!,nullfunc!,nullfunc!),
                0.0, FlatND(), zeros(ndim), zeros(ndim), zeros(ndim), zeros(ndim)
        )
end
