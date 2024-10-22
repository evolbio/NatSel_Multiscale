using MLS, Plots

########################################################
## Create Figure 2 in the manuscript
## Dynamics of zbar over patches ignoring group size, shows w/in group dyn
#  Copy code here and past in Julia REPL or other command line interface
using MLS, Plots
N = 1000;
t_incr = 100;
κ = 100;

qbar = 0.9;
s = 0.9;
r_target = 0.95;
_, x1, M, r = d_qbar(N, qbar, κ, s, r_target);
zbarw, pl = cycle_dynamics(qbar, N, M, x1, κ, s, r, κ, 3; t_incr=t_incr);

#y = zeros(4)
s_array = reverse(push!(collect(0.1:0.2:0.5),0.9));
for i in 1:4
	zbarw, _ = cycle_dynamics(qbar, N, M, x1, κ, s_array[i], r, κ, 3;
					t_incr=t_incr, show_plot=false)
	#y[i] = zbarw[1,2]
	scatter!(pl,[1.005],[zbarw[1,2]],color=mma[i],marker=:circle,markersize=5,
			markerstrokecolor=mma[i])
	annotate!(pl, 0.85, zbarw[1,2], text(string(s_array[i]), 10))
	#plot!(pl, (1:length(zbarw[:]))/(t_incr+1), zbarw[:], color=mma[i], lw=2)
end
display(pl)

# Adjust path for output for local computer
savefig(pl, "/Users/steve/Desktop/temporalDyn.pdf")

########################################################
## Various tests related to dynamics and calculations used for
## making figure

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
