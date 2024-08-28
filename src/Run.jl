using MLS, Plots

q0 = 0.1;
m1 = 0.1;
m2 = 0;
dt = 1.00;
T = 100;

t,q = discrete_dq(q0,m1,m2,dt,T);
plot(t,q)

sol = continuous_dq(q0,m1,m2,T);
plot!(sol)

