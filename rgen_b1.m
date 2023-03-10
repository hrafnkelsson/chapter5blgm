function b1_out=rgen_b1(tau,s2_tau,kappa2_1_tau,phi_1_tau)

global Z2 Z2mZ2 Nx2 Ny2 tau0 M2inv H2

inv_cov_b1 = sparse((Z2mZ2)/s2_tau + (M2inv - phi_1_tau*H2)/kappa2_1_tau);
mean_b1 = inv_cov_b1\((Z2'*(tau-tau0))/s2_tau);
Bchol = chol(inv_cov_b1);
b1_out = mean_b1 + Bchol'\randn(Nx2*Ny2,1);
B0 = ones(1,Nx2*Ny2);
vec0=inv_cov_b1\B0';
b1_out=b1_out-vec0*inv(B0*vec0)*(B0*b1_out);