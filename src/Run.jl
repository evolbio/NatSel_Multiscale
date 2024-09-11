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

#############################

q0 = 0.1;
T = 1;
m1 = 1/T;
r = 0.1;

out = [MLS.qbar_prime(q0,1/T,m1+d,T,r) - q0 for d in -0.1:0.01:0.1];

#############################

q0 = 0.1;
s = 1;
r = 0.1;
k = 1;	# take k as e^t for period t of competition w/in groups, t = log k
m1 = log(k*(1-r)/(1-r+r*s));

out = [MLS.qbar_prime(q0,m1,m1+d,s,r,k) - q0 for d in -0.1:0.01:0.1];

