function pipost = postdens_dim2a(theta)

% The negative log-posterior density based on the exponential distribution
% The sigma parameter varies with longitude.

global Q1 Z1 K y lam1 lam2 Nx1 Ny1

Ktheta = (Nx1 - 2)*(Ny1 - 2);
logs0 = log(0.05);

sigvec = zeros(K,1);
sigvec(1:Nx1) = logs0*ones(1,Nx1);
sigvec((K - Nx1 + 1):K) = logs0*ones(1,Nx1);
sigvec(1 + Nx1*(1:(Ny1 - 2))) = logs0*ones(1,(Ny1 - 2));
sigvec(Nx1 + Nx1*(1:(Ny1  - 2))) = logs0*ones(1,(Ny1 - 2));
sigvec((sigvec==0)) = theta(1:Ktheta);

n = length(y);
pipostvec = zeros(1,n);
%xi = -0.4 + 0.4*exp(theta(Kgx+1))/(1 + exp(theta(Kgx+1)));
%sigvec = Qchol*[exp(theta(1)); (exp(theta(K + 2))*z)];
sig = exp(Z1*sigvec);

for j = 1:n
    pipostvec(j) = explike([sig(j)]',y(j));    
end
pipost = sum(pipostvec) + 0.5*(exp(-2*theta(Ktheta + 1)))*sigvec'*Q1*sigvec + lam2*exp(theta(Ktheta + 1)) - theta(Ktheta + 1); % + 0.5*(1/0.1^2)*(mean(sigvec) - 0.2)^2
end