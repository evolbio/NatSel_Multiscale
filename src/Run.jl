using MLS, Plots

## Global frequency change over groups, shows that simple calculated
#  optimum increases against all variants (full parameter space not checked)

N = 1000;
qbar = 0.2;
κ = 5;
s = 0.2;
r_target = 0.8;
out, x1, M, r = d_qbar(N, qbar, κ, s, r_target; stepsz = 0.01, pts = 2);

out, x1, M, r = d_qbar(N, qbar, κ, s, r_target; stepsz = 0.1, pts = 5);

# dynamics
N = 1000;
qbar = 0.2;
κ = 5;
s = 0.2;
r_target = 0.8;
_, x1, M, r = d_qbar(N, qbar, κ, s, r_target);
dyn = dynamics(qbar, N, M, x1, x1+0.1, s, r, κ, 15);

# polymorphism
dynamics(qbar, N, M, x1-0.09599944802, x1+0.1, s, r, κ, 100000)[end]-qbar
# 8.461650824465039e-11
dynamics(qbar, N, M, x1-0.8749685628, x1+1.1, s, r, κ, 100000)[end]-qbar
# 3.7109093575793395e-11

# start at different qbar and see where things end up
dynamics(0.7, N, M, x1-0.8749685628, x1+1.1, s, r, κ, 100000)[end]
# converges to 0.20000000003718574, suggested that is an attractor
dynamics(0.03, N, M, x1-0.8749685628, x1+1.1, s, r, κ, 100000)[end]
# converges to 0.20000000003714702
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
