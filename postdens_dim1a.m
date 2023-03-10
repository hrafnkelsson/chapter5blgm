function pipost = postdens_dim1a(theta)

% The negative log-posterior density based on the exponential distribution
% The sigma parameter varies with longitude. 

global Rx X Kgx y lam1 lam2 TE

%if (theta(Kgx + 2) <= 0)
%    pipost = Inf
%else
    n = length(y);
    pipostvec = zeros(1,n);
    z = theta(2:Kgx);
    sigvec = TE*[theta(1); (exp(theta(Kgx + 1))*z)];
    %sigvec = TE*[theta(1); ((theta(Kgx + 1))*z)];
    sig = exp(X*sigvec);
    for j = 1:n
        pipostvec(j) = explike([sig(j)]',y(j));    
    end
    %pipost = sum(pipostvec) + 0.5*(z'*z) + lam1*exp(theta(1)) - theta(1) + lam2*(theta(Kgx + 2)) + theta(Kgx+1) + 2*log(1 + exp(-theta(Kgx+1)));
    pipost = sum(pipostvec) + 0.5*(z'*z) + lam1*exp(theta(1)) - theta(1) + lam2*exp(theta(Kgx + 1)) - theta(Kgx + 1);
    %pipost = sum(pipostvec) + 0.5*z'*z + lam1*exp(theta(1)) - theta(1) + 0.5*(1/1^2)*(theta(Kgx + 2) - log(0.01))^2;
%end
end
