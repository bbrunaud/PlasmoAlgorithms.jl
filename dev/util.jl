function get_frac(a)
	if length(a) == 1
		if abs(a%1) < 1e-10 || abs(a%1) > 1-1e-10
			return 0.0
		end
		if a > 0
			return a%1
		else
			return 1+(a%1)
		end
	end

	for i in 1:length(a)
		if abs(a[i]%1)<1e-10 || abs(a[i]%1) > 1-1e-10
			a[i] = 0.0
		else
			if a[i] > 0
				a[i] = a[i]%1
			else
				a[i]= 1+(a[i]%1)
			end
		end 
	end
	return a
end