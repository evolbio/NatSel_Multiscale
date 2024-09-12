using MLS, Plots

## Global frequency change over groups, shows that simple calculated
#  optimum increases against all variants (full parameter space not checked)

N = 1000;
qbar = 0.2;
κ = 5;
s = 0.2;
r_target = 0.8;
out, x1, M, r = d_qbar_prime(N, qbar, κ, s, r_target);

########################################################
## Frequency change analysis within a single patch

q0 = 0.1;
m1 = 0.1;
m2 = 0;
dt = 1.00;
T = 100;

t,q = MLS.discrete_dq(q0,m1,m2,dt,T);
plot(t,q)

sol = MLS.continuous_dq(q0,m1,m2,T);
plot!(sol)



