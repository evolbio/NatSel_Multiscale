module MLS
using DifferentialEquations, Distributions, Printf, Plots
export discrete_dq, continuous_dq, qt, d_qbar, dynamics,
		cycle_dynamics

####################################################################
# colors, based on standard Mathematica palette

mma = [RGB(0.3684,0.50678,0.7098),RGB(0.8807,0.61104,0.14204),
			RGB(0.56018,0.69157,0.19489), RGB(0.92253,0.38563,0.20918)];

####################################################################

## Code for calculating qbar', shows predicted optimum is correct 
# p is probability of group with that has freq q
# test corr: calc_corr(corr_binomial(N, M, qbar),(0:N)/N)

# driver to calculate change in qbar relative to predicted optimum
function d_qbar(N, qbar, κ, s, r_target; stepsz = 0.01, pts = 2)
	M, r = find_M(N, r_target)
	x1 = z_star(κ, r, s);
	rng = stepsz * pts
	out = [qbar_prime(qbar,N,M,x1,x1+d,s,r,κ) - qbar for d in -rng:stepsz:rng]
	return out, x1, M, r
end

calc_corr(p,q) = (sum(p.*q.^2) - sum(p.*q)^2) / (sum(p.*q)*(1-sum(p.*q)))
corr_binomial(N, M, qb) = 
	[qb*pdf(Binomial(N-M,qb),k-M)+(1-qb)*pdf(Binomial(N-M,qb),k) for k in 0:N]
# x_i = e^m_i, assuming t = 1 throughout
y_local(q, x1, x2) = q*x1 + (1 .- q)*x2
z_star(κ, r, s) = κ*(1-r)/(1-r+r*s)

function find_M(N::Int, target_r::Float64)
	approx_M = 0.5 * (1 + sqrt(1 - 4 * N + 4 * N^2 * target_r))
	M = Int(round(approx_M))
	r = (N + M^2 - M) / N^2
	return M, r 
end

function qbar_prime(qbar, N, M, x1, x2, s, r, κ)
	p = corr_binomial(N, M, clamp(qbar,0,1))	# freq of groups with q
	q = collect(0:N) / N						# values of q for groups
	y = y_local(q, x1, x2)
	ybar = sum(p.*y)
	wq = (κ .- y).^s
	wbar = sum(p .* wq)
	α = x1 - x2
	q_prime = @. q + α * q * (1 - q) / y
	return sum(p .* q_prime .* wq) / wbar
end

function dynamics(qbar, N, M, x1, x2, s, r, κ, steps)
	history = zeros(steps+1)
	history[1] = qbar
	for i in 2:(steps+1)
		history[i] = qbar_prime(history[i-1], N, M, x1, x2, s, r, κ)
	end
	return history
end

function cycle_dynamics(qbar, N, M, x1, x2, s, r, κ, steps;
			t_incr=100, show_plot=true)
	# get the frequencies at the start of each cycle
	q0 = dynamics(qbar, N, M, x1, x2, s, r, κ, steps)
	# col of matrix is the ave within grp trait over time for given cycle
	zbarw = zeros(t_incr+1,steps+1)
	tsteps = 0:(1/t_incr):1
	qi = collect(0:N) / N			# values of q for groups
	αt = [x1^t - x2^t for t in tsteps]
	yt = [q*x1^t + (1-q)*x2^t for t in tsteps, q in qi]
	qt = [qi[j] + (αt[i]*qi[j]*(1-qi[j]))/yt[i,j]
			for j in 1:(N+1), i in 1:(t_incr+1)]
	zt = [qt[j,i]*x1 + (1-qt[j,i])*x2
			for j in 1:(N+1), i in 1:(t_incr+1)]
	for j in 1:steps+1
		p = corr_binomial(N, M, clamp(q0[j],0,1))	# freq of groups with qi
		for i in 1:(t_incr+1)
			zbarw[i,j] = sum(p .* zt[:,i])
		end
	end
	pl = nothing
	if show_plot
		pl = plot((1:length(zbarw[:]))/(t_incr+1), zbarw[:], legend=:none,
			color=mma[1], lw=2, xlabel="Rounds of within-group competition",
			ylabel="Average trait value within groups", guidefontsize=13,
			tickfontsize=11)
		display(pl) 
	end
	return zbarw, pl
end

# Calculate freq change in haploid model with two competitors
# Can use for within group changes

s(m1, m2) = m1 - m2
s(m1, m2, dt) = exp(m1*dt) - exp(m2*dt)
wbar(q, m1, m2) = 1 + q*m1 + (1-q)*m2
wbar(q, m1, m2, dt) = q*exp(m1*dt) + (1 - q)*exp(m2*dt)

dq(q, m1, m2, dt) = s(m1, m2, dt) * q*(1-q) / wbar(q, m1, m2, dt)
dq(q, m1, m2) = s(m1, m2) * q * (1-q)
qt(q, m1, m2, t) = discrete_dq(q,m1,m2,t,t)[2][2]

function discrete_dq(q,m1,m2,dt,T)
	steps = 0:dt:T
	qval = zeros(length(steps))
	qval[1] = q
	for i in 2:length(steps)
		q += dq(q,m1,m2,dt)
		qval[i] = q
	end
	return steps, qval
end

function continuous_dq(q0,m1,m2,T)
	prob = ODEProblem((u,p,t) -> dq(u,m1,m2), q0, (0,T))
	sol = solve(prob, Tsit5(), reltol = 1e-12, abstol = 1e-12)
end

end # module MLS
