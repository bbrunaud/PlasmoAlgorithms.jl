function cross(bgraph::Plasmo.PlasmoGraph, lgraph::Plasmo.PlasmoGraph; 
	heuristic=:lagfirst, 
	max_iterations_lag::Int64=60, 
	max_iterations_benders::Int64=100,
	benders_cuts=[:LP],
	ϵ=1e-5,
	benders_UBupdatefrequency=1,
	benders_timelimit=3600,
	benders_verbose=false, 
	benders_LB=NaN, 
	is_nonconvex=false)

end