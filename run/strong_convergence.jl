using Plots			# Plots
using GR

using Distributions # Random Numbers!
using DataFrames    # Read TSV Files!
using StatsBase     # Statistics
using Gallium		# Debugger
# mimic "import wsc"
include("../src/wsc.jl")

mkpath("../results/")

function time_error(eq, dt_axis;err_type="weak_error")
	path = "../results/weak_strong_error/$(eq)/"
	mkpath(path)

	# M = 10^4 for ABM
	# M = 10^5 for GBM
	M = 10^4
	norm = 1/M;

	weak_axis = []
	strong_axis = []

	for dt in dt_axis
		fin_t = 1;
		EM_end = []
		s_dif = []
		X_true = nothing
		X_aprox = nothing
		dBt = nothing
		strong_error = 0

		# time variables
		time_axis,dt,
		frame_time,
		frame_axis,
		frame_index_axis = wsc.get_time_axis(N,0,fin_t,0;dt=dt)

		if eq == "GBM"
			for m in 1:M
				dBt, Bt = wsc.get_dBt(dt, time_axis)
				X_aprox = wsc.EM_GBM_aprox(dBt, Bt, a, b, dt, X_0, time_axis)
				X_true = wsc.GBM_sol(dBt, Bt,a, b, dt, X_0, time_axis)
				#@show X_aprox[1], X_true[1]

				push!(s_dif, abs(X_true[end] - X_aprox[end]))
			end

			strong_error = mean(s_dif)

		elseif eq == "ABM"
			for m in 1:M
				dBt, Bt = wsc.get_dBt(dt, time_axis)
				X_aprox = wsc.EM_ABM_aprox(dBt, Bt, a, b, dt, X_0, time_axis)
				X_true = wsc.ABM_sol(dBt, Bt,a, b, dt, X_0, time_axis)
				#@show X_aprox[1], X_true[1]

				push!(s_dif, abs(X_true[end] - X_aprox[end]))
			end

			strong_error = mean(s_dif)
		end

		push!(strong_axis, strong_error)

		Plots.plot(X_true, label="an_sol")
		Plots.plot!(X_aprox, label="em_apr")
		Plots.png(path*"$(dt).png")
	end
	fit_l, title = wsc.f_line(strong_axis, dt_axis, "linreg")
	wsc.p_err(path,dt_axis,strong_axis,fit_l,err_type,title,true)
end

# length of solution array
N = 2^8

# Eq parameters
a = 0.2
b = 0.01

X_0 = 1

# dt values
dt = 1/(2^11)
s_dt_axis = [dt, 2*dt, 4*dt, 8*dt, 16*dt, 32*dt, 64*dt]


#@time time_error("GBM", s_dt_axis, err_type="strong_error")
@time time_error("ABM", s_dt_axis,err_type="strong_error")
