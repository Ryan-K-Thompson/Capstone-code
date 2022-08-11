syms N %sigma_prime_f epsilon_prime_f sigma_max epsilon_amp b c E 
tic
sigma_prime_f = 1464;
epsilon_prime_f = 0.262;
sigma_max = 384;
epsilon_amp = 0.004;
b = -0.143;
c = -0.619;
E = 71*10^3;

eqn = sigma_max*epsilon_amp == (sigma_prime_f^2)/E*((2*N)^(2*b)) + sigma_prime_f*epsilon_prime_f*((2*N)^(b+c));
% 
% N = solve(eqn, N)

N = vpasolve(eqn,N)
toc