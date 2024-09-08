module MLS
using DifferentialEquations
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

end # module MLS
