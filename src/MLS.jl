module MLS
using DifferentialEquations, Distributions, Integrals, Printf
export discrete_dq, continuous_dq, qt

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

function qbar_prime(qbar, m1, m2, s, r, k)
	g = (1-r)/r
	a = g*qbar
	b = g*(1-qbar)
	T = 1
	prob = IntegralProblem((q,p) -> pdf(Beta(a,b),q)*wbar(q,m1,m2,T), (0,1))
	soln = solve(prob, QuadGKJL())
	ybar = soln.u
	soln.resid > 1e-7 ? println("ybar resid = ", soln.resid) : nothing
	prob = IntegralProblem((q,p) -> 
		pdf(Beta(a,b),q)*discrete_dq(q,m1,m2,T,T)[2][2]*(k-wbar(q,m1,m2,T))^s,
		(0,1))
	soln = solve(prob, QuadGKJL())
	soln.resid > 1e-7 ? println("qbar_prime resid = ", soln.resid) : nothing
	return soln.u / (k-ybar)^s
end

############

function qbar_prime_old(qbar, m1, m2, T, r)
	g = (1-r)/r
	a = g*qbar
	b = g*(1-qbar)
	# set k so that m* = s = 1/T
	k = ℯ*(1-r+r/T)/(1-r)
	prob = IntegralProblem((q,p) -> pdf(Beta(a,b),q)*wbar(q,m1,m2,T), (0,1))
	soln = solve(prob, QuadGKJL())
	ybar = soln.u
	soln.resid > 1e-7 ? println("ybar resid = ", soln.resid) : nothing
	prob = IntegralProblem((q,p) -> 
			pdf(Beta(a,b),q)*discrete_dq(q,m1,m2,T,T)[2][2]*(k-wbar(q,m1,m2,T)),
			(0,1))
	soln = solve(prob, QuadGKJL())
	soln.resid > 1e-7 ? println("qbar_prime resid = ", soln.resid) : nothing
	return soln.u / (k-ybar)
end

end # module MLS
