using Plots

export f_line, p_err

function f_line(error_axis, dt_axis, p_type)
	dt_axis = convert(Array{Float64,1},dt_axis)
	dt_log	 = map(x -> Base.Math.log(x), dt_axis)
	error_axis = convert(Array{Float64,1},error_axis)
	error_log	 = map(x -> Base.Math.log(x), error_axis)

	if p_type == "linreg"
		intercept, slope = linreg(dt_log,error_log)
		f_fit(x) = exp(intercept)*(x^slope)
		line = map(x->f_fit(x),dt_axis)
		title = "fit slope: $(slope)"

		return line, title

	elseif p_type == "lineEq"
		slope = (error_log[end] - error_log[1])/(dt_log[end] - dt_log[1])
		f_fit(x, s, err) = s*(x-dt_log[1]) + err[1]
		line = map(x->f_fit(x,slope,error_log),dt_axis)
		title = "fit slope: $(slope)"

		return line, title
	end
end

function p_err(path,dt_axis,error_axis,fit_l,e_type,title,log)
	Plots.plot(dt_axis,error_axis,
				xlabel="dt",
				ylabel=e_type,
				label=e_type)
	Plots.scatter!(dt_axis,error_axis,
				title=title,
				label="error at dt")
	Plots.plot!(dt_axis,fit_l, label="fit")
	if log
		Plots.plot!(xaxis=:log,yaxis=:log)
	end
	slope = 0
	if e_type == "weak_error"
		slope = 1
	elseif e_type == "strong_error"
		slope = 0.5
	end
	f_decay(x) = x^slope
	decay_line = map(x->f_decay(x),dt_axis)

	Plots.plot!(dt_axis,decay_line,
				label="decay $(slope)",
				line=:dash)

	Plots.png(path*"$(e_type).png")
end
