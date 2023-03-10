% PSHA spatial model in one dimension.

clear all
clear global
close all

global Rx X Kgx y lam1 lam2 TE

% load the data 
[magn,long,lat] = readvars('PanzCat_Mc2_Jun.xlsx','Range','A2:C665');
plot(long,lat,'k.')
indloc = find((lat < 64.15).*(lat > 63.85).*(long < -19.9).*(long > -22.4));
hold, plot(long(indloc),lat(indloc),'r.')
hold

magnb = magn(indloc);
longb = long(indloc);
latb = lat(indloc);

sigm1 = 0.7447;
%y = -sigm1*log(1 - rand(2000,1));
xi = -0.05;
y = magnb - 1.99;
ysort = sort(y);
n = length(y);
[gppar,gpparci] = gpfit(y)
[exppar,expparci] = expfit(y)
sigm = exppar; %mean(y);
p = (1:n)/(n + 1);
t = 0.01:0.01:4.5;
ft = exp(-t/sigm);
ft2 = (1 + gppar(1)*t/gppar(2)).^(-1/gppar(1));
%plot(-log(1 - p),ysort,'k.')
figure(2), plot(-sigm*log(1 - p),ysort,'*')
hold
plot(t,t,'k-')
hold

figure(3), plot(t,1 - ft,'b-')
hold
plot(t,1 - ft2,'r-')
plot(ysort,p,'k*')
hold

% Y < sigma/(-xi)

% (-xi)*ymax = sigma

xiplot = -0.35:0.01:0;
sigmaplot = (-xiplot)*max(y) ;

figure(4),plot(xiplot,sigmaplot,'b-')
hold
plot(gppar(1),gppar(2),'r*')
plot(gpparci(:,1),gppar(2)*ones(2,1),'r-')
plot(gppar(1)*ones(2,1),gpparci(:,2),'r-')
hold

% The spatial dimensions of the rectangular area.
min_lon = -22.5;
max_lon = -19.8;

% The number of equally spaced interior knots.
kx = 6;
% kx = 11;
% kx = 16;
% kx = 21;
% kx = 26;

% Delta x and Delta y.
dx = 1/(kx+1);

% The order of the splines.
M = 4;

% Determine the number of functions.
Nx1 = kx + M;

% The epsilon-knots
epsilon_x = dx*[0:(kx+1)];

% the tau-knots.
tau_x = zeros(1,kx+2*M);
tau_x(1:M) = epsilon_x(1)*ones(1,M);
tau_x(M+1:kx+M) = epsilon_x(2:kx+1);
tau_x(kx+M+1:kx+2*M) = epsilon_x(kx+2)*ones(1,M);

% Vector with values of x and y.
x_bsp = (longb - min_lon)/(max_lon - min_lon);
lx = length(x_bsp);

% Compute the x-splines and the y-splines.
X = spline_functions(x_bsp,tau_x,dx,kx,M);

% Compute R1
Kgx = Nx1;
D1 = diag([1 2*ones(1,Nx1-2) 1]);
Rx = zeros(1,Kgx);
Rx(2) = -1;
Rx = toeplitz(Rx) + D1;

% Evaluate the likelihood

%Dens = @(t)-DensEval_gplrc_var1_c1_constraints(t,RC);

TE = Pgen_new_eps(Kgx);
lam1 = 1/0.5;
lam2 = 1/(0.1/Kgx);
% postdens_dim1([log(0.48); randn(Kgx - 1,1); 0.5; log(0.001)]) %  Z
postdens_dim1b([log(0.2); randn(Kgx - 1,1); 0.47; log(0.001)])   % non Z

% Z
[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim1b(theta),[log(0.2); randn(Kgx - 1,1); 0.47; log(0.001)]);
% non Z
%[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim1(theta),[log(0.48); randn(Kgx - 1,1); 0.5; log(0.001)]);
postdens_dim1b(t_m)


sigvec = TE*[exp(t_m(1)); (exp(t_m(Kgx + 2))*t_m(2:Kgx))]
sig = X*sigvec
xi = t_m(Kgx+1);
sig_nu = exp(t_m(Kgx + 2))
for j = 1:n
    pipostvec(j) = gplike([xi, sig(j)]',y(j));    
end
loglikelihood = -sum(pipostvec)

spar = sqrt(diag(inv(H)));
spar(1) 
spar(Kgx + 1)
spar(Kgx + 2)

[sort_long,indx_long] = sort(longb);
figure(1), plot(longb(indx_long),sig(indx_long),'k-')
xlabel('Longitude')
ylabel('Sigma')
title('Sigma (black)')
print(1, '-dpdf', 'sigma_oef_dim1.pdf')

%figure(5), plot(sort(y./(sig/(-xi))),'k*')
%figure(6), plot(sort(longb),5+(sig(indx_long)/(-xi)),'k-')
%hold
figure(6), plot(sort(longb),2+(sig(indx_long)*(2^(xi)-1)/(xi)),'r-')
hold
plot(sort(longb),2 + y(indx_long),'b*')
hold
xlabel('Longitude')
ylabel('Magnitude')
title('Magnitude data (blue), median (red)')
print(6, '-dpdf', 'magnitude_vs_longitude_oef_dim1.pdf')

mean(y(indx_long(1:165)))
mean(y(indx_long(166:330))) 
mean(y(indx_long(331:495)))    
mean(y(indx_long(496:663)))     

%mean(y(ind_long(1:17)))/max(y(ind_long(1:17))) = 0.4894 
%mean(y(ind_long(18:34)))/max(y(ind_long(18:34))) = 0.3684
%mean(y(ind_long(35:52)))/max(y(ind_long(35:52))) = 0.3686
%mean(y(ind_long(53:68)))/max(y(ind_long(53:68))) = 0.3479
