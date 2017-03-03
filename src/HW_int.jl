
module HW_int

	const A_SOL = 4  # analytic solution


	# question 1 b) 

	using FastGaussQuadrature
	using Roots
	using Sobol
	using Plots
	using Distributions
	using NullableArrays
	using DataStructures  # OrderedDict

	# choose plots backend
	plotlyjs()

	# demand function
	q(p) = 2*(p.^-0.5)

	# gauss-legendre adjustment factors for map change
	ba2(lb,ub) = (ub-lb)/2
	ab2(lb,ub) = (lb+ub)/2

	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply, 
	# i.e. domestic + export demand = 2
	function dd(p,t1,t2)
	    exp(t1)/p .+ exp(t2) .* p^-0.5 - 2
	end

	function fzero_wrap(f,lb,ub) 
		try
			fzero(f,lb,ub)
		catch
			println("no eqm price found")
			return(NaN)
		end
	end

	print_f(x)  = abs(100*(A_SOL-x))/A_SOL
	print_fn(x) = round(print_f(x),5)

	function plot_q1()

		# run for all n
		d = OrderedDict()
		ns = [10;100;1000]
		for n in ns
			d[n] = Dict()
			d[n][:laguerre] = question_1b(n)
			d[n][:MC] = question_1c(n)
			d[n][:QMC] = question_1d(n)
		end

		# build 2x3 plot matrix
		# 			Laguerre 	MC 		QMC					
		# points    pts vs vals at n = 10
		# error     n vs error
		l = @layout grid(2,3)
		# pts locations for n=10
		p1 = Any[]
		for (k,v) in d[10]
			p = scatter(v[:x],v[:y],legend=false,xlims=(0.5,4.5),ylims=(0.5,2.5),title=k)
			if k==:laguerre
				yaxis!("Eval Points")
			end
			push!(p1,p)
		end
		# errors vs n
		for (k,v) in d[10]
			p = scatter(ns,[print_f(d[ns[1]][k][:I]);print_f(d[ns[2]][k][:I]);print_f(d[ns[3]][k][:I])],xaxis=(:log10),legend=false,xlims=(9,2*ns[3]),xticks=ns,xformatter=x->round(Int,x),yformatter=:scientific,m=(2,:cross,:red))
			if k==:laguerre
				yaxis!("Percent Error")
			end
			push!(p1,p)
		end
		plot(p1...,layout=l)
		savefig("HW-int.pdf")

	end


	function question_1b(n)

		gl = gausslegendre(n);

		# bounds on integration
		a = 1
		b = 4

		# use transformation formula to map into [-1,1]
		pts  = ba2(a,b).*gl[1] .+ ab2(a,b)	# integration points
		vals = q(pts)  # function values at those points
 		Integ = ba2(a,b) * (gl[2]' * vals ) # do integration

		# plot
		# plots only for n=10
		# if n==10
		# 	figure()
		# 	plot(pts,vals,"o-")
		# 	title("Gauss Laguerre")
		# end

		println("estimated change in CS using $n gauss legendre nodes is $(Integ)")
		println("i.e. an error of $(print_fn(Integ)) percent")
		println("")
		return Dict(:y=>vals, :x=>pts, :I => Integ)
	end

	# question 1 c)

	function question_1c(n)

		# get n random numbers from [1,4]
		pts = rand(n)*3 + 1	
		vals = q(pts)

		# integrate: Monte carlo is defined for the "hypercube" [0,1]
		# we need to adjust the "volume" of this cube to be 3
		# Integ = 3*mean(vals) - 2
		Integ = 3*mean(vals) 
		
		# plot
		# if n==10
		# 	figure()
		# 	plot(pts,vals,"o")
		# 	ylim([0;5])
		# 	axhline(mean(vals),color="red")
		# 	title("Monte Carlo")
		# end

		println("estimated change in CS using $n monte carlo nodes is $Integ)")
		println("i.e. an error of $(print_fn(Integ)) percent")
		println("")
		return Dict(:y=>vals, :x=>pts, :I => Integ)

	end

	function question_1d(n)

		s = SobolSeq(1,[1],[4])  # 1-dimensional sobol sequence in [1,pstar]
		pts = zeros(n)
		for i in 1:n
			pts[i] = next(s)[1]
		end

		vals = q(pts)

		# integrate: Monte carlo is defined for the hypercube [0,1]
		# we need to extend the length of that interval to be 3
		Integ = 3*mean(vals)
		
		# plot
		# if n==10
		# 	figure()
		# 	plot(pts,vals,"o")
		# 	axhline(mean(vals),color="red")
		# 	title("Quasi Monte Carlo")
		# end

		println("estimated change in CS using $n Quasi monte carlo nodes is $Integ)")
		println("i.e. an error of $(print_fn(Integ)) percent")
		println("")
		return Dict(:y=>vals, :x=>pts, :I => Integ)

	end

	# question 2

	function question_2a(n)

		gh = gausshermite(n)

		Sigma = hcat([0.02, 0.01],[0.01,0.01])
		Omega = chol(Sigma)

		mu = [0.0;0.0]

		# kronecker product of grids and weights
		gr = hcat(kron(ones(n),gh[1]),kron(gh[1],ones(n)))
		wt = kron(gh[2],gh[2]) / pi	# watch out for the pi!

		# make adjustment for correlation in shocks
		grids = Omega * gr'	 + zeros(2,n*n)   # zeros here would be a matrix with mu

		# find eqm price at each combination of theta1,theta2
		pstar = NullableArray(Float64,n*n)  # create an array of nullable for float64
		for i in 1:length(pstar)
			try
				pstar[i] = fzero(x->dd(x,grids[1,i],grids[2,i]),1.0)
			catch
				# assign nothing to pstar[i]:
				# no eqm price found
				# hence this is a missing value
			end
		end

		# plot
		if n==10
			p = scatter(grids[1,:],grids[2,:])
			title!("Question 2a: Gauss hermite theta grid")
		end
		wt = wt[!pstar.isnull]
		pstar = dropnull(pstar)
		EP = dot(wt,pstar)  # using all non-null values to compute this
		VAR = dot(wt, (pstar .- EP).^2 )

		if n==10
			return Dict("E[p]"=>EP, "Var[p]"=>VAR, "plot"=>p)
		else
			return Dict("E[p]"=>EP, "Var[p]"=>VAR)

		end

	end

	function question_2b(n)

		# for fairness, let's also create n^2 points as in 2a
		n = n^2

		Sigma = hcat([0.02, 0.01],[0.01,0.01])
		M = MvNormal(Sigma)	# create mean zero joint normal distribution
		pts = rand(M,n)	# just draw from it randomly

		# find eqm price at each combination of theta1,theta2
		pstar = NullableArray(Float64,n*n)  # create an array of nullable for float64
		for i in 1:length(pstar)
			try
				pstar[i] = fzero(x->dd(x,grids[1,i],grids[2,i]),1.0)
			catch
				# assign nothing to pstar[i]:
				# no eqm price found
				# hence this is a missing value
			end
		end
		
		# plot
		if n==100
			p = scatter(pts[1,:],pts[2,:])
			title!("Question 2b: Monte Carlo theta grid")
		end

		pstar = dropnull(pstar)
		EP = mean(pstar)

		VAR = mean( (pstar .- EP).^2 )
		if n==100
			return Dict("E[p]"=>EP, "Var[p]"=>VAR, "plot"=>p)
		else
			return Dict("E[p]"=>EP, "Var[p]"=>VAR)
		end

	end

	function question_2bonus(n)

		# for fairness, let's also create n^2 points as in 2a
		n = n^2

		s = SobolSeq(2,[1,1],[4,4])  # 2-dimensional sobol sequence 
		pts = hcat([next(s) for i=1:n])
		if n == 100
			println("here's the sobol sequence")
			println(pts)
		end

		# find eqm price at each combination of theta1,theta2
		pstar = NullableArray(Float64,n*n)  # create an array of nullable for float64
		for i in 1:length(pstar)
			try
				pstar[i] = fzero(x->dd(x,grids[1,i],grids[2,i]),0.001,15)
			catch
				# assign nothing to pstar[i]:
				# no eqm price found
				# hence this is a missing value
			end
		end
		println("it's very hard to find the eqm price for those random values")
		
		# plot
		if n==100
			p = scatter([pts[i][1] for i in 1:n],[pts[i][2] for i in 1:n])
			title!("Question 2bonus: Quasi Monte Carlo theta grid")
		end

		pstar = dropnull(pstar)
		EP = mean(pstar)

		VAR = mean( (pstar .- EP).^2 )
		if n==100
			return Dict("E[p]"=>EP, "Var[p]"=>VAR, "plot"=>p)
		else
			return Dict("E[p]"=>EP, "Var[p]"=>VAR)
		end


	end

	# function to run all questions
	function runall()
		info("Running all of HW-integration")
		info("question 1:")
		plot_q1()
		for n in (10,100,1000)
			info("============================")
			info("now showing results for n=$n")
			# info("question 1b:")
			# question_1b(n)	# make sure your function prints some kind of result!
			# info("question 1c:")
			# question_1c(n)
			# info("question 1d:")
			# question_1d(n)
			# println("")
			info("question 2a:")
			if n < 1000
				q2 = question_2a(n)
				println(q2)
			end
			if n == 10
				info("question 2b:")
				q2b = question_2b(n)
				println(q2b)
				info("bonus question: Quasi monte carlo:")
				q2bo = question_2bonus(n)
				println(q2bo)
			else
				info("skipping question 2b: takes too long")
			end

			println()
		end
	end
	info("end of HW-integration")

end

