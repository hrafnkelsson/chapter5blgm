function a1_out=rgen_a1(mu,s2_mu,kappa2_1_mu,phi_1_mu)

global Z1 Z1mZ1 Nx1 Ny1 mu0 M1inv H1

inv_cov_a1 = sparse((Z1mZ1)/s2_mu+(M1inv-phi_1_mu*H1)/kappa2_1_mu);
mean_a1 = inv_cov_a1\((Z1'*(mu-mu0))/s2_mu);
Achol = chol(inv_cov_a1);
a1_out=mean_a1+Achol'\randn(Nx1*Ny1,1);
A0=ones(1,Nx1*Ny1);
vec0=inv_cov_a1\A0';
a1_out=a1_out-vec0*inv(A0*vec0)*(A0*a1_out);
