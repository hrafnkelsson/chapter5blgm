function pipost = postdens_dim1_test(theta)

% The negative log-posterior density based on the GP distribution
% The sigma parameter varies with longitude, xi is unknown

global Rx X Kgx y lam1 lam2 TE xiref x_bsp 

%xiref = 0.75;
%if (theta(Kgx + 2) <= 0)
%    pipost = Inf
%else
    n = length(y);
    pipostvec = zeros(1,n);
    xi = -xiref + xiref*exp(theta(5))/(1 + exp(theta(5)));   
    sig = exp(theta(1) + theta(2)*(x_bsp - 0.5) + theta(3)*(x_bsp - 0.5).^2 + theta(4)*(x_bsp - 0.5).^3);
    for j = 1:n
        j
        pipostvec(j) = gplike([xi, sig(j)]',y(j));    
    end
    %pipost = sum(pipostvec) + 0.5*(z'*z) + lam1*exp(theta(1)) - theta(1) + lam2*(theta(Kgx + 2)) + theta(Kgx+1) + 2*log(1 + exp(-theta(Kgx+1)));
    pipost = sum(pipostvec) + theta(5) + 2*log(1 + exp(-theta(5)));
    %pipost = sum(pipostvec) + 0.5*z'*z + lam1*exp(theta(1)) - theta(1) + 0.5*(1/1^2)*(theta(Kgx + 2) - log(0.01))^2;
%end
end
