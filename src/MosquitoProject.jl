module MosquitoProject

import Base

using CSV, DelimitedFiles, Dates
using DataFrames
using LinearAlgebra
using Interpolations
using JLD, MAT
using Statistics
using DifferentialEquations


export  MyTuple,
        FlatND,
        PbcND,
        Ψ1,
        Ψ2,
        term_mag,
        Polynomial,
        LaguerreFunction,
        LegendreP,
        ChebyshevT,
        ChebyshevU,
        export_clean_traces,
        save_data,
        mat_load_data,
        jld_load_data,
        nullfunc!,
        ODEModels,
        MosquitoModel,
        MosquitoInference,
        init_mosquitoinference,
        sparse_bayesian_fit,
        SBL,
        SBLvar,
        define_sde_prob,
        sde_sim,
        get_rhs!,
        get_w!,
        mynorm

function mynorm(vec)
    vv = 0.0
    for i in eachindex(vec)
        vv += vec[i]*vec[i]
    end
    vv = sqrt(vv)
    return vv
end

mutable struct MyTuple{T1, T2, T3, T4}
    velpot::T1
    visual::T2
    odor::T3
    flow::T4
    MyTuple{T1,T2,T3,T4}(x1, x2, x3, x4) where {T1,T2,T3,T4} = new( convert(T1,x1), convert(T2,x2), convert(T3,x3), convert(T4,x4) ) 
end
MyTuple(x1::T1, x2::T2, x3::T3, x4::T4) where {T1, T2, T3, T4} = MyTuple{T1, T2, T3, T4}(x1, x2, x3, x4)


Base.getindex(D::MyTuple{T1, T2, T3, T4}, sym::Symbol) where {T1, T2, T3, T4} = getfield(D, sym)
Base.setindex!(D::MyTuple{T1, T2, T3, T4}, val, sym::Symbol) where {T1, T2, T3, T4} = try 
                                                                                          setfield!(D, sym, val )
                                                                                      catch any_error 
                                                                                          setfield!(D, sym, convert(typeof(D[sym]),val) )
                                                                                      end
Base.convert(::Type{MyTuple{T1, T2, T3, T4}}, D::MyTuple) where {T1, T2, T3, T4} = MyTuple{T1,T2,T3,T4}(
    convert(T1, D.velpot), convert(T2, D.visual), convert(T3, D.odor), convert(T4, D.flow)
)
Base.length(::MyTuple{T1, T2, T3, T4}) where {T1, T2, T3, T4} = 4

include("base.jl")
include("clean.jl")
include("load.jl")
include("model.jl")
include("optimization.jl")
include("resim.jl")



end