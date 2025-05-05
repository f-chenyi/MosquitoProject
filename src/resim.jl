function define_sde_prob(model::M; u0=zeros(2*model.ndim), tspan=(0.0, 1.0)) where M <: ODEModels
    prob = SDEProblem(get_rhs!, get_w!, u0, tspan, model);
    return prob
end

function sde_sim(model::M, u0::Vector{<:Number}, tspan::Tuple{<:Number, <:Number}=(0.0, 1.0), solver=EM();
                 dts=0.01, ode_args...) where M <: ODEModels
    prob = define_sde_prob(model; u0=u0, tspan=tspan);
    sol = solve(prob, solver; saveat=dts, tstops=collect(tspan[1]:dts:tspan[2]), ode_args...)
    cache_sol = Array(sol');
    tmax = floor( Int, (prob.tspan[2] - prob.tspan[1]) / dts ) + 1

    return sol.t[1:tmax], cache_sol[1:tmax,1:model.ndim], cache_sol[1:tmax,(model.ndim+1):(2*model.ndim)]
end


function sde_sim(model::M, u0::Vector{Vector{T1}}, tspan::Vector{Tuple{T2, T3}}, solver=EM();
                      dts=0.01, ode_args...) where {M <: ODEModels, T1 <: Number, T2 <: Number, T3 <: Number}
    
    (length(u0) != length(tspan)) && (error("DimensionError: dimensions of u0 and tspan don't match"))
    
    prob = define_sde_prob(model);
    
    rarray = Vector{Matrix{Float64}}(undef, length(u0))
    varray = Vector{Matrix{Float64}}(undef, length(u0))
    tarray = Vector{Vector{Float64}}(undef, length(u0))
    
    for n = 1:length(u0)
        prob_n = remake( prob; tspan=tspan[n], u0=u0[n]);
        sol_n = solve(prob_n, solver; saveat=dts, tstops=collect(prob_n.tspan[1]:dts:prob_n.tspan[2]), ode_args...)
        cache_sol = Array(sol_n');
        tmax = floor( Int, (sol_n.t[end] - sol_n.t[1]) / dts ) + 1
        rarray[n] = cache_sol[1:tmax,1:model.ndim];
        varray[n] = cache_sol[1:tmax,(model.ndim+1):(2*model.ndim)];
        tarray[n] = sol_n.t[1:tmax]
    end

    return tarray, rarray, varray
end






# define ODE equations
function get_rhs!(du, u, t, 
                  ndim::Int, vel::Vector{T}, pos::Vector{T}, bvec::Vector{T}, uvec::Vector{T},
                  islearned::MyTuple{Bool, Bool, Bool, Bool},
                  funcs, coeffs, params, bfields) where T
    
    fill!(du, 0.0)
    fill!(pos, 0.0)
    fill!(vel, 0.0)
    

    # dx/dt = v
    for k = 1:ndim
        du[k]=u[ndim + k]
    end
    
    # vel <- v
    for k=1:ndim
        pos[k] = u[k]
        vel[k] = u[ndim+k]
    end

    # a field
    vmag = mynorm(vel)
    vel ./= vmag

    for sym in fieldnames(typeof(islearned))
        
        (!islearned[sym]) && (continue)
        
        # compute b 
        fill!(bvec, 0.0)
        bmag = bfields[sym](bvec, pos)
        
        # compute a dot b
        ab_dot = 0.0;
        for k in 1:ndim
            ab_dot += vel[k] * bvec[k];
        end
        
        fmag = 0.0
        for idx in eachindex(funcs[sym])
            # fmag
            (typeof(funcs[sym][idx])<:Ψ1) && ( fmag = term_mag(funcs[sym][idx], vmag, params[sym][1]))
            (typeof(funcs[sym][idx])<:Ψ2) && ( fmag = term_mag(funcs[sym][idx], vmag, bmag, ab_dot, params[sym][1], params[sym][2]) ) 
                
            # uvec
            (funcs[sym][idx].uvec == :ahat) && (copy!(uvec, vel))
            (funcs[sym][idx].uvec == :a) && (copy!(uvec, vel); uvec .*= vmag;)
            (funcs[sym][idx].uvec == :bhat) && (copy!(uvec, bvec))
            (funcs[sym][idx].uvec == :b) && (copy!(uvec, bvec); uvec .*= bmag;)
            (funcs[sym][idx].uvec == :bhatorth) && ( uvec .= (bvec .- vel .* ab_dot) ./ vmag )
            (funcs[sym][idx].uvec == :bhatorth2) && ( uvec .= (bvec .- vel .* ab_dot) )
            for k in 1:ndim
                du[ndim+k] +=  coeffs[sym][idx] * fmag * uvec[k]
            end
        end
        
    end
    
    nothing
end

get_rhs!(du, u, M::MosquitoModel, t) = get_rhs!(du, u, t, 
                                       M.ndim, M.velocity, M.position, M.bvec, M.uvec,
                                       M.islearned, M.funcs, M.coeffs, M.params, M.bfields);


function get_w!(dw, u, t, ndim::Int, sqrtD::T) where T
    fill!(dw,0.0)
    for k=1:ndim
        dw[ndim+k] = sqrtD
    end
    nothing
end


get_w!(dw, u, M::MosquitoModel, t) = get_w!(dw, u, t, M.ndim, sqrt.(M.diffD))