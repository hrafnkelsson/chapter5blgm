function pipost = postdens_dim2(theta)

% The negative log-posterior density based on the GP distribution
% The sigma parameter varies with longitude, xi is unknown

global Q1 Z1 K y lam1 lam2

n = length(y);
pipostvec = zeros(1,n);
%xi = -0.4 + 0.4*exp(theta(Kgx+1))/(1 + exp(theta(Kgx+1)));
xi = theta(K+1);
sigvec = theta(1:K);
%sigvec = Qchol*[exp(theta(1)); (exp(theta(K + 2))*z)];
sig = Z1*sigvec;

for j = 1:n
    pipostvec(j) = gplike([xi, sig(j)]',y(j));    
end
pipost = sum(pipostvec) + 0.5*(exp(-2*theta(K + 2)))*sigvec'*Q1*sigvec + 0.5*(1/0.1^2)*(mean(sigvec) - 0.2)^2 + lam2*exp(theta(K + 2)) - theta(K + 2);
end