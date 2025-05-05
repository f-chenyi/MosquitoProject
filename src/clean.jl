parse_date(d::AbstractString) = DateTime(d, dateformat"yyyy-mm-dd HH:MM:SS.ssss")

function export_clean_traces(file;
        min_track_length = 1.0, min_start_time = 300.0, max_start_time=Inf,
        wall_buffer_x = 0.0 , wall_buffer_y = 0.0, wall_buffer_z = 0.0
)

    data = CSV.read(file*".csv", DataFrame)
    df = convert(DataFrame,data)
    groups = groupby(df, 1)
    
    starting_times, clean_traces = save_data(groups, file; 
        min_track_length=min_track_length, min_start_time=min_start_time, max_start_time=max_start_time,
        wall_buffer_x=wall_buffer_x, wall_buffer_y=wall_buffer_y, wall_buffer_z=wall_buffer_z
    )
    
    return starting_times, clean_traces
    
end

function save_data(groups, file; 
        min_track_length = 1.0, min_start_time = 300.0, max_start_time = Inf,
        wall_buffer_x = 0.0 , wall_buffer_y = 0.0, wall_buffer_z = 0.0
    )
    # first, compute positions of wall
    min_x = Vector{Float64}()
    max_x = Vector{Float64}()
    min_y = Vector{Float64}()
    max_y = Vector{Float64}()
    min_z = Vector{Float64}()
    max_z = Vector{Float64}()
    for g in groups
        push!(min_x,minimum(g[:,3]))
        push!(max_x,maximum(g[:,3]))
        push!(min_y,minimum(g[:,4]))
        push!(max_y,maximum(g[:,4]))
        push!(min_z,minimum(g[:,5]))
        push!(max_z,maximum(g[:,5]))
    end

    wall_robust = [(quantile(min_x,0.01),quantile(max_x,0.99)),(quantile(min_y,0.01),quantile(max_y,0.99)),(quantile(min_z,0.01),quantile(max_z,0.99))]
    println("Wall x position: ", wall_robust[1])
    println("Wall y position: ", wall_robust[2])
    println("Wall z position: ", wall_robust[3])

    # Filter tracks based on time requirements
    experiment_start_time = parse_date(groups[1].datetime[1][1:end-3])
    long_traces = []
    for g in groups
        track_end_time = parse_date(g.datetime[end][1:end-3])
        track_start_time = parse_date(g.datetime[1][1:end-3])
        if Dates.value(track_end_time - track_start_time)/1000 > min_track_length && Dates.value(track_start_time-experiment_start_time)/1000 > min_start_time && Dates.value(track_start_time-experiment_start_time)/1000 < max_start_time
            push!(long_traces,g)
        end
    end

    # Filter tracks based on position from wall
    final_traces = []
    times = []
    starting_times = []
    for g in long_traces
        if all(g[:,3] .> wall_robust[1][1]+wall_buffer_x) && all(g[:,3] .< wall_robust[1][2]-wall_buffer_x) && all(g[:,4] .> wall_robust[2][1]+wall_buffer_y) && all(g[:,4] .< wall_robust[2][2]-wall_buffer_y) && all(g[:,5] .> wall_robust[3][1]+wall_buffer_z) && all(g[:,5] .< wall_robust[3][2]-wall_buffer_z)
            push!(final_traces,g)
            tmp = []
            for i in 1:length(g.datetime)
                push!(tmp,Dates.value(parse_date(g.datetime[i][1:end-3])-parse_date(g.datetime[1][1:end-3]))/1000)
            end
            push!(starting_times,Dates.value(parse_date(g.datetime[1][1:end-3])-experiment_start_time)/1000)
            push!(times,tmp)
        end
    end


    clean_traces = Vector{Matrix{Float64}}()
    for g in 1:length(final_traces)
        if sum(diff(times[g]) .<= 0.0) == 0
            ts = LinRange(0.0,times[g][end],Int(round(times[g][end]*100)+1))
            interp_linear_xs = linear_interpolation(times[g], final_traces[g][:,3])
            interp_linear_ys = linear_interpolation(times[g], final_traces[g][:,4])
            interp_linear_zs = linear_interpolation(times[g], final_traces[g][:,5])

            temp = zeros(length(ts),5)
            temp[:,1] = interp_linear_xs.(ts) / 100 # convert cm to m
            temp[:,2] = interp_linear_ys.(ts) / 100
            temp[:,3] = interp_linear_zs.(ts) / 100
            temp[:,4] = ts

            for j in 1:length(ts)
                if round(ts[j],digits=2) in times[g]
                    temp[j,5] = 1.0
                else
                    temp[j,5] = 0.0
                end
            end
            push!(clean_traces,temp)
        else
            println("Skipping track ",g)
            deleteat!(starting_times, g)
        end
    end

#     for i in 1:length(clean_traces)
#         clean_traces[i][:,4] .+= starting_times[i]
#     end

    save(file*"_clean_traces.jld", "clean_traces", clean_traces)
    return starting_times, clean_traces

end