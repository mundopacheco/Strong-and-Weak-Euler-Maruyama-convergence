export EM_GBM_aprox, GBM_sol,
ABM_sol, get_dBt, get_time_axis


"""
Function that computes the Euler-Maruyama aproximation
to the Geometric Brownian Motion equation
"""
function EM_GBM_aprox(dBt, Bt,a, b, dt, X_0, time_axis)
		X_aprox = []
		X_i = X_0
		push!(X_aprox,X_i)
    for t in eachindex(time_axis[2:end])
        X_i += (a*X_i)*dt + (b*X_i)*dBt[t]
		push!(X_aprox,X_i)
    end
    return X_aprox
end

"""
Function that computes the solution
to the Geometric Brownian Motion equation
"""
function GBM_sol(dBt, Bt, a, b, dt, X_0, time_axis)
		X_sol = []
    for t in eachindex(time_axis)
        X_t = X_0 * exp((a-((1/2)*(b^2)))*time_axis[t] + (b*Bt[t]))
		push!(X_sol,X_t)
    end
    return X_sol
end

"""
Function that computes the Euler-Maruyama aproximation
to the Aritmetic Brownian Motion equation
"""
function EM_ABM_aprox(dBt, Bt,a, b, dt, X_0, time_axis)
		X_aprox = []
		X_i = X_0
		push!(X_aprox,X_i)
    for t in eachindex(time_axis[2:end])
        X_i += (a)*dt + (b)*dBt[t]
		push!(X_aprox,X_i)
    end
    return X_aprox
end

"""
Function that computes the solution
to the Aritmetic Brownian Motion equation
"""
function ABM_sol(dBt, Bt, a, b, dt, X_0, time_axis)
		X_sol = []
    for t in eachindex(time_axis)
        X_t = X_0 + a*time_axis[t] + b*Bt[t]
		push!(X_sol,X_t)
    end
    return X_sol
end


function get_dBt(dt, time_axis)
	dBt = randn(length(time_axis))*sqrt(dt)
	Bt = cumsum(dBt)
	#Plots.plot(dBt, label="dBt")
	#Plots.plot!(Bt, label="Bt")
	#Plots.png("../results/weak_strong_error/noise.png")
	return dBt, Bt
end

function get_time_axis(N,init_t,fin_t,frame_amount;dt=0)

	if dt == 0
		dt = 1/(N)
	end
	# time variables

	time_axis	 		= init_t:dt:fin_t
	frame_time   		= length(time_axis)/frame_amount
	frame_axis   		= []
	frame_index_axis	= []

	#@show frame_time,frame_amount
	t = init_t
	netxt_frame = 1

	for i in 1:length(time_axis)
		if i >= netxt_frame
			push!(frame_index_axis, i)
			netxt_frame += frame_time
		end
	end

	frame_index_axis[end] = length(time_axis)
	#@show length(time_axis)
	if frame_amount != 0
		for i in 1:frame_amount
			push!(frame_axis,time_axis[frame_index_axis[i]])
		end
	end

	return time_axis,dt,frame_time,
		frame_axis,frame_index_axis
end
