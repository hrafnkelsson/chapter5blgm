function pipost = postdens_dim0(theta)

% The negative log-posterior density based on the GP distribution
% The sigma parameter is unknown but constant with longitude, xi is unknown

global y lam1 xiref alp bet

n = length(y);
pipostvec = zeros(1,n);
xi = -xiref + xiref*exp(theta(2))/(1 + exp(theta(2)));
sigma = exp(theta(1));
for j = 1:n
    pipostvec(j) = gplike([xi, sigma]',y(j));   
end
%pipost = sum(pipostvec) + 0.5*(z'*z) + lam1*exp(theta(1)) - theta(1) + lam2*(theta(Kgx + 2)) + theta(Kgx+1) + 2*log(1 + exp(-theta(Kgx+1)));
pipost = sum(pipostvec) + lam1*exp(theta(1)) - theta(1) + theta(2) + 2*log(1 + exp(-theta(2))) - (alp - 1)*log(xi/xiref + 1) - (bet - 1)*log(-xi/xiref); 
%pipost = sum(pipostvec) + 0.5*z'*z + lam1*exp(theta(1)) - theta(1) + 0.5*(1/1^2)*(theta(Kgx + 2) - log(0.01))^2;
end
