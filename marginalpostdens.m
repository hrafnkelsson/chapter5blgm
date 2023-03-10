function pitheta = marginalpostdens(theta)

global Q1 Z1 Qdata muvec lambda M

Qcond = (1/theta^2)*Q1 + Z1'*Qdata*Z1;
mucond = Qcond\(Z1'*Qdata*muvec);
L = chol(Qcond);
logdet = sum(log(diag(L)));
muQmu = mucond'*Qcond*mucond; 
pitheta = - lambda*theta - (M - 1)*log(theta) - logdet + 0.5*muQmu;
end

