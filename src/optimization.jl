abstract type ModelInferencePipeline end
mutable struct MosquitoInference <: ModelInferencePipeline
    pos_track::Vector{Matrix{Float64}}  # x_i
    vel_track::Vector{Matrix{Float64}}  # v_i
    tobs::Vector{Array{Float64}}       # t_i
    model::MosquitoModel          # model with functional term
    Θ::Vector{Array{Float64,3}}   # F(x, v)
    dVdt::Vector{Matrix{Float64}} # dVdt - F_learned 
    coef_full::Vector{Float64}    # w_μ
    bigind::BitVector
    G::Matrix{Float64}
    b::Vector{Float64}            # dV/dt - F_learned
    cacheG::Matrix{Float64}
    cacheb::Vector{Float64}
    
    MosquitoInference(robs::Vector{Matrix{Float64}}, vobs::Vector{Matrix{Float64}}, 
                      tobs::Vector{Array{Float64}}, model::MosquitoModel) = new(
        robs, vobs, tobs, model,
        Vector{Array{Float64,3}}(undef, length(robs)),
        Vector{Array{Float64,3}}(undef, length(robs)),
        [0.0], trues(1), zeros(1,1), zeros(1), zeros(1,1), zeros(1)
    )
end



# ============= FUNCTIONS for the inference pipeline ==================

# ------------- Initialization -------------
function init_mosquitoinference(P::MosquitoInference)

    # get keys 
    allkeys, learnkey = internal_dim_check(P)
    nterm = sum( [length(P.model.funcs[learnkey_]) for learnkey_ in learnkey] )
    
    # init Θ (using the first slice to store the velocity trajectories)
    # init dVdt
    for i in eachindex(P.Θ)
        P.Θ[i] = zeros(size(P.vel_track[i])...,nterm)
        # dts = P.tobs[i][2]-P.tobs[i][1]
        dts = @view(P.tobs[i][2:end,:]) .- @view(P.tobs[i][1:end-1,:]) 
        P.dVdt[i] = (@view(P.vel_track[i][2:end,:]) .- @view(P.vel_track[i][1:end-1,:]) ) ./ dts
    end
    
    # update dVdt
    learnedkeys = allkeys[[P.model.islearned[key] for key in allkeys]]
    for lkey in learnedkeys
        println("learned force: $(lkey)")
        _update_dVdt!(P.dVdt, P.pos_track, P.vel_track, P.model.funcs[lkey], P.model.params[lkey], P.model.coeffs[lkey], 
                      lkey, P.model.ndim, P.model.bfields[lkey])
    end
    
    # init coeffs
    P.coef_full = zeros(nterm)
    
    # init bigind
    P.bigind = trues(nterm)
      
    # init b
    _init_b!(P)
    
    # init G
    P.G = zeros(length(P.b), nterm)
    
    return nothing
end



function internal_dim_check(P::MosquitoInference)

    # check the input data has the same dimension has the model
    any(size.(P.pos_track,2) .!= P.model.ndim) && (error("DimensionError: the dimension of the model must be the same as that of the positional data"));
    any(size.(P.vel_track,2) .!= P.model.ndim) && (error("DimensionError: the dimension of the model must be the same as that of the velocity data"));
    (length(P.pos_track) != length(P.vel_track)) && (error("DimensionError: pos and vel data must have the same number of particles"))
    
    
    # check only one force term in in the learning process
    keys = [fieldnames(typeof(P.model.islearning))...]
    islearning_vec = [P.model.islearning[k] for k in keys]
    (sum(islearning_vec) != 1) && (error("need at least one and only one force term to be inferred")) 
    
    return keys, keys[islearning_vec]
    
end




function _update_dVdt!(dVdt::Vector{Matrix{Float64}}, pos::Vector{Matrix{Float64}}, vel::Vector{Matrix{Float64}},
                       funcs, params::Vector{Float64}, coeffs::Vector{Float64}, key::Symbol, ndim::Int, get_bfield!::Function)
    
    any(size.(dVdt,1) .!= (size.(pos,1).-1) ) && (error("DimensionError: size(dVdt,1) == size(pos,1)-1 not satisfied"))
    
    ri_vec = zeros(ndim)
    vi_vec = zeros(ndim)
    u_vec  = zeros(ndim)
    b_vec  = zeros(ndim)
    
    for pid in eachindex(pos)
        
        for t in 1:size(dVdt[pid],1)

            for k in 1:ndim
                ri_vec[k] = pos[pid][t,k]
                vi_vec[k] = vel[pid][t,k]
            end
            
            # a_mag
            vi_mag = mynorm(vi_vec)
            vi_vec ./= vi_mag
            
            # b_mag
            b_mag = get_bfield!(b_vec, ri_vec);
            ab_dot = 0.0;
            for k in 1:ndim
                ab_dot += vi_vec[k] * b_vec[k];
            end
#             ab_dot = dot(vi_vec, b_vec);
            
            fmag = 0.0            
            for idx in eachindex(funcs)
                (typeof(funcs[idx])<:Ψ1) && (fmag = term_mag(funcs[idx], vi_mag, params[1]))
                (typeof(funcs[idx])<:Ψ2) && ( fmag = term_mag(funcs[idx], vi_mag, b_mag, ab_dot, params[1], params[2]) ) 
                (funcs[idx].uvec == :ahat) && (copy!(u_vec, vi_vec))
                (funcs[idx].uvec == :a) && (copy!(u_vec, vi_vec); u_vec .*= vi_mag;)
                (funcs[idx].uvec == :bhat) && (copy!(u_vec, b_vec))
                (funcs[idx].uvec == :b) && (copy!(u_vec, b_vec); u_vec .*= b_mag;)
                (funcs[idx].uvec == :bhatorth) && ( u_vec .= (b_vec .- vi_vec .* ab_dot) ./ vi_mag )
                (funcs[idx].uvec == :bhatorth2) && ( u_vec .= (b_vec .- vi_vec .* ab_dot) )
                for k in 1:ndim
                    dVdt[pid][t,k] -= coeffs[idx] * fmag * u_vec[k]
                end
            end
            
        end

    end
    
    return nothing
end


function _init_b!(P::MosquitoInference)
    row_array = length.(P.dVdt)
    P.b = zeros(sum(row_array))
    count = 0 
    for ip in eachindex(P.dVdt)
        copyto!(@view(P.b[count+1:count+row_array[ip]]),
                P.dVdt[ip])
        count += row_array[ip]
    end
    return nothing
end



function build_theta!(Theta_s::Vector{Array{Float64,3}}, pos::Vector{Matrix{Float64}}, vel::Vector{Matrix{Float64}},
                       funcs,  params::Vector{Float64}, key::Symbol, ndim::Int, get_bfield!::Function)
    
    ri_vec = zeros(ndim)
    vi_vec = zeros(ndim)
    u_vec  = zeros(ndim)
    b_vec  = zeros(ndim)

    for pid in eachindex(pos)

        fill!(Theta_s[pid], 0.0)
        
        for t in 1:size(pos[pid],1)

            for k in 1:ndim
                ri_vec[k] = pos[pid][t,k]
                vi_vec[k] = vel[pid][t,k]
            end
            
            # a_mag
            vi_mag = mynorm(vi_vec)
            vi_vec ./= vi_mag
            
            # b_mag
            b_mag = get_bfield!(b_vec, ri_vec);
            ab_dot = 0.0;
            for k in 1:ndim
                ab_dot += vi_vec[k] * b_vec[k];
            end
#             ab_dot = dot(vi_vec, b_vec);
            
            fmag = 0.0            
            for idx in eachindex(funcs)
                # fmag
                (typeof(funcs[idx])<:Ψ1) && ( fmag = term_mag(funcs[idx], vi_mag, params[1]))
                (typeof(funcs[idx])<:Ψ2) && ( fmag = term_mag(funcs[idx], vi_mag, b_mag, ab_dot, params[1], params[2]) ) 
                
                # uvec
                (funcs[idx].uvec == :ahat) && (copy!(u_vec, vi_vec))
                (funcs[idx].uvec == :a) && (copy!(u_vec, vi_vec); u_vec .*= vi_mag;)
                (funcs[idx].uvec == :bhat) && (copy!(u_vec, b_vec))
                (funcs[idx].uvec == :b) && (copy!(u_vec, b_vec); u_vec .*= b_mag;)
                (funcs[idx].uvec == :bhatorth) && ( u_vec .= (b_vec .- vi_vec .* ab_dot) ./ vi_mag )
                (funcs[idx].uvec == :bhatorth2) && ( u_vec .= (b_vec .- vi_vec .* ab_dot) )
                for k in 1:ndim
                    Theta_s[pid][t, k, idx] += fmag * u_vec[k]
                end
            end
            
        end

    end
end


@inline _update_theta!(P::MosquitoInference, model::MosquitoModel, sym::Symbol) = build_theta!(P.Θ, P.pos_track, 
                       P.vel_track, model.funcs[sym], model.params[sym], sym, model.ndim, model.bfields[sym])



function _update_G!(G::Matrix{Float64}, Θ::Vector{Array{Float64,3}})
    
    count=0

    for ip in eachindex(Θ)
        
        L1 = size(Θ[ip],1) - 1
        L2 = size(Θ[ip],2)
        nrows= L1*L2
        
        for n in axes(G,2)
            copyto!( @view( G[count+1:count+nrows,n] ),
                     @view( Θ[ip][1:end-1,:,n])
                )
        end
        
        count += nrows

    end

end



function sparse_bayesian_fit(self::MosquitoInference, params_new, sym::Symbol; opt_args...)
    
    copy!(self.model.params[sym], params_new)
    _update_theta!(self, self.model, sym)
    _update_G!(self.G, self.Θ)
    sbl_res = SBL(self.G, self.b; opt_args...)
    
    return sbl_res

end



function SBL(Phi::Matrix{Float64}, Y::Vector{Float64}; 
            MAX_ITERS=100, EPSILON=1e-6, lambda=1.0, gamma=0.5*ones(size(Phi,2)))

    # fitting Y = Phi * mu + GaussianNoise(lambda)
    
    N, M = size(Phi)
    @assert N == length(Y)
    
    mu = zeros(M)
    mu_old = zeros(M)
    
    PhiT_Phi = Phi'*Phi;
    G_inv = Diagonal(zeros(M));
    Sigma = zeros(M,M)
    Xi = zeros(M, N);
    deltaY = similar(Y)
    
    count = 0
    
    while true
        
        copyto!(mu_old, mu)
        
        G_inv.diag .= 1 ./ gamma
        
        Sigma .= PhiT_Phi ./ lambda .+ G_inv .+ 1e-8
        
        Q = lu!(Sigma);
        Sigma .= Q \ I;
        
#         LinearAlgebra.inv!(cholesky!(Sigma))
#         Sigma .= inv(Sigma)
        
        
        BLAS.gemm!('N','T',1/lambda, Sigma, Phi, false, Xi)
        mul!( mu, Xi, Y)
        
        copyto!(deltaY, Y)
        BLAS.gemm!('N','N',-1.0, Phi, mu, true, deltaY)
        
        for k in eachindex(gamma)
            gamma[k] = mu[k]*mu[k] + Sigma[k,k]
        end
        
        lambda = sum(abs2, deltaY)
        # lambda /= (N-sum(gamma))
        lambda /= N
                
        count += 1
        (count >= MAX_ITERS) && (break)
        (all(abs.(mu_old .- mu) .< EPSILON)) && (break)
        
    end
    
    return  1/2*log(lambda), mu

end



function SBLvar(Phi::Matrix{Float64}, Y::Vector{Float64}; 
            MAX_ITERS=100, EPSILON=1e-6, lambda=1.0, gamma=0.5*ones(size(Phi,2)))

    # fitting Y = Phi * mu + GaussianNoise(lambda)
    
    N, M = size(Phi)
    @assert N == length(Y)
    
    mu = zeros(M)
    mu_old = zeros(M)
    
    PhiT_Phi = Phi'*Phi;
    G_inv = Diagonal(zeros(M));
    Sigma = zeros(M,M)
    Xi = zeros(M, N);
    deltaY = similar(Y)
    
    count = 0
    
    while true
        
        copyto!(mu_old, mu)
        
        G_inv.diag .= 1 ./ gamma
        
        Sigma .= PhiT_Phi ./ lambda .+ G_inv .+ 1e-8
        
        Q = lu!(Sigma);
        Sigma .= Q \ I;
        
#         LinearAlgebra.inv!(cholesky!(Sigma))
#         Sigma .= inv(Sigma)
        
        
        BLAS.gemm!('N','T',1/lambda, Sigma, Phi, false, Xi)
        mul!( mu, Xi, Y)
        
        copyto!(deltaY, Y)
        BLAS.gemm!('N','N',-1.0, Phi, mu, true, deltaY)
        
        for k in eachindex(gamma)
            gamma[k] = mu[k]*mu[k] + Sigma[k,k]
        end
        
        lambda = sum(abs2, deltaY)
        # lambda /= (N-sum(gamma))
        lambda /= N
                
        count += 1
        (count >= MAX_ITERS) && (break)
        (all(abs.(mu_old .- mu) .< EPSILON)) && (break)
        
    end
    
    return  1/2*log(lambda), mu, Sigma

end