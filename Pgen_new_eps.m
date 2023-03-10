function P = Pgen_new_eps(Kgx)
P = zeros(Kgx,Kgx);
for j = 1:(Kgx)
   P(j,1:j) = ones(1,j);   
end
end