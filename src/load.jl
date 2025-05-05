function mat_load_data(; data_dir = "./", file_name = "MTracks_trap_on_fan_off_odor_off_crop_100.mat")
    
    matdata = matread(data_dir * file_name)

    _pos = convert(Vector{Matrix{Float64}}, matdata["pos"][:,1])
    _vel = convert(Vector{Matrix{Float64}}, matdata["vel"][:,1])
    _t = convert(Vector{Array{Float64}}, matdata["t"][:,1]);

    return _pos, _vel, _t
end


function jld_load_data(; data_dir = "./", file_name = "tmp.jld",
                         trace_key = "clean_traces", diff_mode = :forward, 
                         len_scale = 1.0, time_scale=1.0, interp = false, dts=0.01, ndim=3)
    
    jlddata = load(data_dir * file_name)
    trdata = jlddata[trace_key]
    
    __t = [trdata[i][:, ndim+1] for i in eachindex(trdata)]
    __pos = [trdata[i][:, 1:ndim] for i in eachindex(trdata) ]

    # discard any trace that is unsorted
    i_sort = findall( issorted.(__t) );
    __t = __t[i_sort]
    __pos = __pos[i_sort]
    
    if interp
        for i in eachindex(__t)
            tn = collect(__t[i][1]:dts:__t[i][end])
            yn = interp1( __t[i], __pos[i], tn)
            __t[i] = tn
            __pos[i] = yn
        end
    end
    
    if diff_mode == :forward
        _t = [__t[i][1:end-1] * time_scale for i in eachindex(__t)]
        _pos = [__pos[i][1:end-1,:] * len_scale for i in eachindex(__pos)]
        _vel = [ diff(__pos[i], dims=1) ./ diff(__t[i], dims=1) * (len_scale/time_scale) for i in eachindex(__t)]
    elseif diff_mode == :central
        _t = [__t[i][2:end-1] * time_scale for i in eachindex(__t)]
        _pos = [__pos[i][2:end-1,:] * len_scale for i in eachindex(__pos)]
        _vel = [ (__pos[i][3:end,:] .- __pos[i][1:end-2,:]) ./ (__t[i][3:end] .- __t[i][1:end-2] ) * (len_scale/time_scale) for i in eachindex(__t)]
    elseif diff_mode == :central2
        _t = [__t[i][3:end-2] * time_scale for i in eachindex(__t)]
        _pos = [__pos[i][3:end-2,:] * len_scale for i in eachindex(__pos)]
        _vel = [ (__pos[i][5:end,:] .- __pos[i][1:end-4,:]) ./ (__t[i][5:end] .- __t[i][1:end-4] ) * (len_scale/time_scale) for i in eachindex(__t)]
    elseif diff_mode == :central3
        _t = [__t[i][4:end-3] * time_scale for i in eachindex(__t)]
        _pos = [__pos[i][4:end-3,:] * len_scale for i in eachindex(__pos)]
        _vel = [ (__pos[i][7:end,:] .- __pos[i][1:end-6,:]) ./ (__t[i][7:end] .- __t[i][1:end-6] ) * (len_scale/time_scale) for i in eachindex(__t)]
    else
        error("diff_mode $(diff_mode) is invalid")
    end

    _t = convert(Vector{Array{Float64}}, _t);
    
    return _pos, _vel, _t
end



function interp1(xpt::AbstractVector{T}, ypt::AbstractMatrix{T}, x::AbstractVector{T}) where {T}

    y = zeros(length(x), size(ypt,2));
    
    for k in axes(y, 2)
        intf = interpolate((xpt,), ypt[:,k], Gridded(Linear()))
        copy!( @view(y[:,k]), intf.(x) );
    end
    
    return y
end
