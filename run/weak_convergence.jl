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
		weak_error = 0

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

				push!(EM_end, X_aprox[end])
			end

			weak_error = abs(X_0*exp(a*fin_t) - norm*sum(EM_end))

		elseif eq == "ABM"
			for m in 1:M
				dBt, Bt = wsc.get_dBt(dt, time_axis)
				X_aprox = wsc.EM_ABM_aprox(dBt, Bt, a, b, dt, X_0, time_axis)
				X_true = wsc.ABM_sol(dBt, Bt,a, b, dt, X_0, time_axis)
				#@show X_aprox[1], X_true[1]

				push!(EM_end, X_aprox[end])
			end

			weak_error = abs( ((a*fin_t)+X_0) - norm*sum(EM_end) )
		end

		push!(weak_axis, weak_error)

		Plots.plot(X_true, label="an_sol")
		Plots.plot!(X_aprox, label="em_apr")
		Plots.png(path*"$(dt).png")
	end
	fit_l, title = wsc.f_line(weak_axis, dt_axis, "linreg")
	wsc.p_err(path,dt_axis,weak_axis,fit_l,err_type,title,true)
end

# length of solution array
N = 2^8

# Eq parameters
a = 0.2
b = 0.01

X_0 = 1

# dt values
w_dt_axis = [1/(2^10),1/(2^9),1/(2^8),1/(2^7),1/(2^6)]


@time time_error("ABM", w_dt_axis, err_type="weak_error")			# "ABM", "GBM"
#time_error("GBM", w_dt_axis, err_type="weak_error")
