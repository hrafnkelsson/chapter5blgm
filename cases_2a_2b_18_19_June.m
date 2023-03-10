% PSHA spatial model in one dimension.

clear all
clear global
close all

global Q1 Z1 K y lam1 lam2 Nx1 Ny1

% 'Sheet1'
% '7-17Jun'   ???
% '18-19Jun'   exponential mmin = 2
% '20-21Jun'   exponential mmin = 2
% '22-23Jun'   exponential mmin = 2
% '24-30Jun'

% load the data 
[magn,long,lat] = readvars('PanzCat_Mc2_Jun_version_2.xlsx','Sheet','18-19Jun','Range','A2:C285');
%[magn,long,lat] = readtable("PanzCat_Mc2_Jun_version_2.xlsx","Sheet","7-17Jun","Range","A2:C285");
%[magn,long,lat] = xlsread('PanzCat_Mc2_Jun_version_2.xlsx',6,'A2:C285');
plot(long,lat,'k.')

mmin = 1.999;
indloc = find((lat < 64.2).*(lat > 63.7).*(magn > mmin));
plot(long(indloc),lat(indloc),'r.')
magnb = magn(indloc);
longb = long(indloc);
latb = lat(indloc);

y = magnb - mmin;
ysort = sort(y);
n = length(y);

% The spatial dimensions of the rectangular area.
min_lon = -22.5 - (22.5-20)/26;
max_lon = -20.0 + (22.5-20)/26;

min_lat = 63.85 - (64.15 - 63.85)/7;
max_lat = 64.15 + (64.15 - 63.85)/7;

% The number of equally spaced interior knots.
kx = 25;    % no. x B-splines = 29
ky = 8;    % no. y B-splines = 12

% Delta x and Delta y.
dx = 1/(kx+1);

% Delta x and Delta y.
dy = 1/(ky+1);

% The order of the splines.
M = 4;

% Determine the number of functions.
Nx1 = kx + M;
Ny1 = ky + M;

% The epsilon-knots
epsilon_x = dx*[0:(kx+1)];
epsilon_y = dy*[0:(ky+1)];

% the tau-knots.
tau_x = zeros(1,kx+2*M);
tau_x(1:M) = epsilon_x(1)*ones(1,M);
tau_x(M+1:kx+M) = epsilon_x(2:kx+1);
tau_x(kx+M+1:kx+2*M) = epsilon_x(kx+2)*ones(1,M);

tau_y = zeros(1,ky+2*M);
tau_y(1:M) = epsilon_y(1)*ones(1,M);
tau_y(M+1:ky+M) = epsilon_y(2:ky+1);
tau_y(ky+M+1:ky+2*M) = epsilon_y(ky+2)*ones(1,M);

% Vector with values of x and y.
x_bsp = (longb - min_lon)/(max_lon - min_lon);
y_bsp = (latb - min_lat)/(max_lat - min_lat);
lx = length(x_bsp);
ly = length(y_bsp);

% Compute the x-splines and the y-splines.
X = spline_functions(x_bsp,tau_x,dx,kx,M);
Y = spline_functions(y_bsp,tau_y,dy,ky,M);

% Create the matrix Z1 with the tensor product of X and Y.
N = n;
K = Nx1*Ny1;
Z1 = zeros(N,Nx1*Ny1);
for i = 1:Ny1
    for j = 1:Nx1
        Z1(:,Nx1*(i-1)+j) = X(:,j).*Y(:,i);
    end
end
Z1 = sparse(Z1);

%-------------------------------------------------------------------------

% Z2 : a matrix with high-dimensional B-spline tensor-products
% Z2 can used to draw figure of the field.

% The number of equally spaced interior knots.
% kx = 47;
% ky = 22;

% Delta x and Delta y.
dx = 1/(kx+1);
dy = 1/(ky+1);

% The order of the splines.
M = 4;

% Determine the number of functions.
Nx2 = kx + M;
Ny2 = ky + M;

% The epsilon-knots
epsilon_x = dx*[0:(kx+1)];
epsilon_y = dy*[0:(ky+1)];

% the tau-knots.
tau_x = zeros(1,kx+2*M);
tau_x(1:M) = epsilon_x(1)*ones(1,M);
tau_x(M+1:kx+M) = epsilon_x(2:kx+1);
tau_x(kx+M+1:kx+2*M) = epsilon_x(kx+2)*ones(1,M);
tau_y = zeros(1,ky+2*M);
tau_y(1:M) = epsilon_y(1)*ones(1,M);
tau_y(M+1:ky+M) = epsilon_y(2:ky+1);
tau_y(ky+M+1:ky+2*M) = epsilon_y(ky+2)*ones(1,M);

% Vector with values of x and y.
sx_temp = min_lon + (max_lon - min_lon)*(0:(1/100):1);
sy_temp = min_lat + (max_lat - min_lat)*(0:(1/32):1);
lx_temp = length(sx_temp);
ly_temp = length(sy_temp);
N2 = lx_temp*ly_temp;
x_bsp_temp = kron(ones(ly_temp,1),sx_temp');
y_bsp_temp = kron(sy_temp',ones(lx_temp,1));

x_bsp = (x_bsp_temp - min_lon)/(max_lon - min_lon);
y_bsp = (y_bsp_temp - min_lat)/(max_lat - min_lat);
N2 = length(x_bsp);

% Compute the x-splines and the y-splines.
X = spline_functions(x_bsp,tau_x,dx,kx,M);
Y = spline_functions(y_bsp,tau_y,dy,ky,M);

% Create the matrix Z2 with the tensor product of X and Y.
Z2 = zeros(N2,Nx2*Ny2);
for i = 1:Ny2
    for j = 1:Nx2
        Z2(:,Nx2*(i-1)+j) = X(:,j).*Y(:,i);
    end
end
Z2 = sparse(Z2);

% Define the matrix Q1.
% H1, C1 and M1
Kgx = Nx1;
Kgy = Ny1;
D1 = diag([1 2*ones(1,Nx1-2) 1]);
D2 = diag([1 2*ones(1,Ny1-2) 1]);
Rx = zeros(1,Kgx);
Rx(2) = -1;
Rx = toeplitz(Rx) + D1;
Ry = zeros(1,Kgy);
Ry(2) = -1;
Ry = toeplitz(Ry) + D2;
H1 = kron(eye(Kgy),Rx)+kron(Ry,eye(Kgx));
Q1 = sparse(H1);
% Qchol = chol(Q1);
% diffQ=Q1 - Qchol'*Qchol;
% sum(sum(diffQ))

% Evaluate the likelihood  Case 2b

%Dens = @(t)-DensEval_gplrc_var1_c1_constraints(t,RC);

TE = Pgen_new_eps(Kgx);
lam1 = 1/0.5;
lam2 = 1/(1/Kgx);
Ktheta = (Nx1 - 2)*(Ny1 - 2);
% postdens_dim1([log(0.48); randn(Kgx - 1,1); -0.5; log(0.001)]) %  Z
% postdens_dim1([log(0.48); randn(Kgx - 1,1); 0.5; log(0.001)])   % non Z
postdens_dim2a([(log(0.2)+0.01*randn(Ktheta,1)); log(0.001)])   % non Z

stop

%postdens_dim1([log(0.48); randn(Kgx - 1,1); 0.9; log(0.001)])
% Z
[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim2a(theta),[(log(0.196)+0.01*randn(Ktheta,1)); log(0.001)]);
% non Z
%[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim1(theta),[log(0.48); randn(Kgx - 1,1); 0.5; log(0.001)]);
postdens_dim2a(t_m)

logs0 = log(0.05);
sigvec = zeros(K,1);
sigvec(1:Nx1) = logs0*ones(1,Nx1);
sigvec((K - Nx1 + 1):K) = logs0*ones(1,Nx1);
sigvec(1 + Nx1*(1:(Ny1 - 2))) = logs0*ones(1,(Ny1 - 2));
sigvec(Nx1 + Nx1*(1:(Ny1  - 2))) = logs0*ones(1,(Ny1 - 2));
sigvec((sigvec==0)) = t_m(1:Ktheta);
sig = exp(Z1*sigvec);
sigplot = exp(Z2*sigvec);

Stheta = H\eye(length(t_m));
corrmat = diag(sqrt(1./diag(Stheta)))*Stheta*diag(sqrt(1./diag(Stheta)));
stderr = sqrt(1./diag(Stheta))

%sigvec = TE*[t_m(1); (exp(t_m(Kgx + 2))*t_m(2:Kgx))]
sigvec = TE*[t_m(1); (exp(t_m(Kgx + 2))*t_m(2:Kgx))]
sig = exp(X*sigvec)
xi = -0.6 + 0.6*exp(t_m(Kgx+1))/(1 + exp(t_m(Kgx+1)))
sig_nu = exp(t_m(Kgx + 2))
for j = 1:n
    pipostvec(j) = gplike([xi, sig(j)]',y(j));    
end
loglikelihood = -sum(pipostvec)

spar = sqrt(diag(inv(H)));
spar(1) 
spar(Kgx + 1)
spar(Kgx + 2)

% Metropolis step for the unknown parmeters
K = length(t_m);
uc = 2.38*sqrt(1/K);
L = 25000;
burnin = 2000;
Lb = burnin + L;
m = 4;  
mut = t_m;
R = chol(H);
theta01 = mut + 0.1*uc*(R\randn(K,1));
theta02 = mut + 0.1*uc*(R\randn(K,1));
theta03 = mut + 0.1*uc*(R\randn(K,1));
theta04 = mut + 0.1*uc*(R\randn(K,1));
Stheta = H\eye(length(t_m));
theta1 = zeros(Lb,K,m);
theta1(1,:,1) = theta01;
theta1(1,:,2) = theta02;
theta1(1,:,3) = theta03;
theta1(1,:,4) = theta04;
for ms = 1:m
for j = 2:Lb
    [ms,j]
    thetastar = theta1(j - 1,:,ms) + uc*(R\randn(K,1))';           % Step 1 
    logpoststar = postdens_dim1(thetastar');
    if logpoststar == Inf 
         theta1(j,:,ms) = theta1(j - 1,:,ms);
    else
        logr = - logpoststar + postdens_dim1((theta1(j - 1,:,ms))');    % Step 2
        if log(rand(1)) < logr
            theta1(j,:,ms) = thetastar;  
        else 
            theta1(j,:,ms) = theta1(j - 1,:,ms);
        end % Step 3.
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figures

minlong = -23;
maxlong = -19;
longplot = 0:0.005:1;
longplot = minlong + longplot*(maxlong - minlong);

% Vector with values of x and y.
x_bsp_plot = (longplot - min_lon)/(max_lon - min_lon);
lx_plot = length(x_bsp_plot);

% Compute the x-splines and the y-splines.
Xplot = spline_functions(x_bsp_plot,tau_x,dx,kx,M);

mu = 5;
[sort_long,indx_long] = sort(longb);

sigma1 = zeros(Lb,lx_plot,m);
for ms = 1:m
    for j = 1:Lb
        [ms,j]
        z = (theta1(j,2:Kgx,ms))';
        sigvec = TE*[theta1(j,1,ms); (exp(theta1(j,Kgx + 2,ms))*z)];
        sig = exp(Xplot*sigvec);
        sigma1(j,:,ms) = sig;  
    end
end

sigmahat = zeros(1,lx_plot);
for k = 1:lx_plot
sigma2a = sigma1((burnin + 1):Lb,k,1);
sigma2b = sigma1((burnin + 1):Lb,k,2);
sigma2c = sigma1((burnin + 1):Lb,k,3);
sigma2d = sigma1((burnin + 1):Lb,k,4);
sigma2mcmc = [sigma2a; sigma2b; sigma2c; sigma2d];
sigmahat(k) = median(sigma2mcmc);
end

pnr = Kgx + 1;
xiref = 0.75;
xi2a = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,1))./(1 + exp(theta1((burnin + 1):Lb,pnr,1)));
xi2b = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,2))./(1 + exp(theta1((burnin + 1):Lb,pnr,2)));
xi2c = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,3))./(1 + exp(theta1((burnin + 1):Lb,pnr,3)));
xi2d = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,4))./(1 + exp(theta1((burnin + 1):Lb,pnr,4)));
xi2mcmc = [xi2a; xi2b; xi2c; xi2d];
xihat = median(xi2mcmc)

%figure(4), hist(xi2mcmc,100)
%figure(1), plot(longplot,sigmahat,'k-')

% Generalized Pareto distribution, spatially varying sigma
minlong = -23;
maxlong = -19;
mu = 5;
[sort_long,indx_long] = sort(longb);

figure(1), plot(longplot,(mu+(sigmahat./(-xihat))),'k-')
hold
plot(longplot,(mu+(sigmahat*((1 - 0.05)^(-xihat) - 1)/(xihat))),'k--')
plot(longplot,(mu+(sigmahat*((1 - 0.25)^(-xihat) - 1)/(xihat))),'k--')
plot(longplot,(mu+(sigmahat*((1 - 0.5)^(-xihat) - 1)/(xihat))),'k--')
plot(longplot,(mu+(sigmahat*((1 - 0.75)^(-xihat) - 1)/(xihat))),'k--')
plot(longplot,(mu+(sigmahat*((1 - 0.95)^(-xihat) - 1)/(xihat))),'k--')
plot(sort(longb),5 + y(indx_long),'k.','MarkerSize',10)
hold
set(gcf,'PaperUnits','centimeters')
set(gcf,'papersize',[15,10])
set(gcf,'paperposition',[0,0,15,10])
set(gcf,'units','centimeters')
set(gcf,'Position',[0,0,15,10])
xlabel('Longitude','interpreter','Latex','fontsize',17)
ylabel('Magnitude','interpreter','Latex','fontsize',17)
set(gca,'fontsize',13)
axs=axis
axis([minlong maxlong mu (mu + 3.0)])
hold off
box on
print(1, '-dpdf', 'gp_spatial_15_magnitude_vs_longitude.pdf')

% Sigma 
sigma0 = zeros(Lb,lx,m);
for ms = 1:m
    for j = 1:Lb
        [ms,j]
        z = (theta1(j,2:Kgx,ms))';
        sigvec = TE*[theta1(j,1,ms); (exp(theta1(j,Kgx + 2,ms))*z)];
        sig = exp(X*sigvec);
        sigma0(j,:,ms) = sig;  
    end
end

sigma0hat = zeros(1,lx);
sigma0mat = zeros(4*L,lx);
for k = 1:lx
sigma0a = sigma0((burnin + 1):Lb,k,1);
sigma0b = sigma0((burnin + 1):Lb,k,2);
sigma0c = sigma0((burnin + 1):Lb,k,3);
sigma0d = sigma0((burnin + 1):Lb,k,4);
sigma0mcmc = [sigma0a; sigma0b; sigma0c; sigma0d];
sigma0mat(:,k) = sigma0mcmc;  
sigma0hat(k) = median(sigma0mcmc);
end


% DIC and pd
Dhat4 = 0;
for i = 1:n
    Dhat4 = Dhat4 + 2*gplike([xihat, sigma0hat(i)],y(i))
end
Dall4 = zeros(1,L*4);
for j = 1:(L*4)
    Dtemp = 0;
    for i = 1:n
        Dtemp = Dtemp + 2*gplike([xi2mcmc(j), sigma0mat(j,i)],y(i));
    end   
    Dall4(j) = Dtemp;
end
Dave4 = mean(Dall4);
pd4 = Dave4 - Dhat4;
DIC4 = Dhat4 + 2*pd4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate the likelihood  Case 2a

%Dens = @(t)-DensEval_gplrc_var1_c1_constraints(t,RC);

TE = Pgen_new_eps(Kgx);
lam1 = 1/0.5;
lam2 = 1/(1/Kgx);
% postdens_dim1([log(0.48); randn(Kgx - 1,1); -0.5; log(0.001)]) %  Z
% postdens_dim1([log(0.48); randn(Kgx - 1,1); 0.5; log(0.001)])   % non Z
postdens_dim1([log(0.48); randn(Kgx - 1,1); 0.9; log(0.001)])   % non Z

%postdens_dim1([log(0.48); randn(Kgx - 1,1); 0.9; log(0.001)])
% Z
[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim1a(theta),[log(0.48); randn(Kgx - 1,1); log(0.001)]);
% non Z
%[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim1(theta),[log(0.48); randn(Kgx - 1,1); log(0.001)]);
postdens_dim1a(t_m)

Stheta = H\eye(length(t_m));
corrmat = diag(sqrt(1./diag(Stheta)))*Stheta*diag(sqrt(1./diag(Stheta)));
stderr = sqrt(1./diag(Stheta))

sigvec = TE*[t_m(1); (exp(t_m(Kgx + 1))*t_m(2:Kgx))]
sig = exp(X*sigvec)
sig_nu = exp(t_m(Kgx + 1))
for j = 1:n
    pipostvec(j) = explike([sig(j)]',y(j));    
end
loglikelihood = -sum(pipostvec)

spar = sqrt(diag(inv(H)));
spar(1) 
spar(Kgx + 1)

% Metropolis step for the unknown parmeters
K = length(t_m);
uc = 2.38*sqrt(1/K);
L = 25000;
burnin = 2000;
Lb = burnin + L;
m = 4;  
mut = t_m;
R = chol(H);
theta01 = mut + 0.1*uc*(R\randn(K,1));
theta02 = mut + 0.1*uc*(R\randn(K,1));
theta03 = mut + 0.1*uc*(R\randn(K,1));
theta04 = mut + 0.1*uc*(R\randn(K,1));
Stheta = H\eye(length(t_m));
theta1a = zeros(Lb,K,m);
theta1a(1,:,1) = theta01;
theta1a(1,:,2) = theta02;
theta1a(1,:,3) = theta03;
theta1a(1,:,4) = theta04;
for ms = 1:m
for j = 2:Lb
    [ms,j]
    thetastar = theta1a(j - 1,:,ms) + uc*(R\randn(K,1))';           % Step 1 
    logpoststar = postdens_dim1a(thetastar');
    if logpoststar == Inf 
         theta1a(j,:,ms) = theta1(j - 1,:,ms);
    else
        logr = - logpoststar + postdens_dim1a((theta1a(j - 1,:,ms))');    % Step 2
        if log(rand(1)) < logr
            theta1a(j,:,ms) = thetastar;  
        else 
            theta1a(j,:,ms) = theta1a(j - 1,:,ms);
        end % Step 3.
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figures

minlong = -23;
maxlong = -19;
longplot = 0:0.005:1;
longplot = minlong + longplot*(maxlong - minlong);

% Vector with values of x and y.
x_bsp_plot = (longplot - min_lon)/(max_lon - min_lon);
lx_plot = length(x_bsp_plot);

% Compute the x-splines and the y-splines.
Xplot = spline_functions(x_bsp_plot,tau_x,dx,kx,M);

mu = 5;
[sort_long,indx_long] = sort(longb);

sigma1a = zeros(Lb,lx_plot,m);
for ms = 1:m
    for j = 1:Lb
        [ms,j]
        z = (theta1a(j,2:Kgx,ms))';
        sigvec = TE*[theta1a(j,1,ms); (exp(theta1a(j,Kgx + 1,ms))*z)];
        sig = exp(Xplot*sigvec);
        sigma1a(j,:,ms) = sig;  
    end
end

sigmahata = zeros(1,lx_plot);
for k = 1:lx_plot
sigma2a = sigma1a((burnin + 1):Lb,k,1);
sigma2b = sigma1a((burnin + 1):Lb,k,2);
sigma2c = sigma1a((burnin + 1):Lb,k,3);
sigma2d = sigma1a((burnin + 1):Lb,k,4);
sigma2mcmc = [sigma2a; sigma2b; sigma2c; sigma2d];
sigmahata(k) = median(sigma2mcmc);
end

% Generalized Pareto distribution, spatially varying sigma
minlong = -23;
maxlong = -19;
mu = 5;
[sort_long,indx_long] = sort(longb);
figure(2), plot(longplot,(mu+(sigmahat*(-log(1 - 0.05)))),'k--'), hold
plot(longplot,(mu+(sigmahat*(-log(1 - 0.25)))),'k--'), 
plot(longplot,(mu+(sigmahat*(-log(1 - 0.5)))),'k--'), 
plot(longplot,(mu+(sigmahat*(-log(1 - 0.75)))),'k--'), 
plot(longplot,(mu+(sigmahat*(-log(1 - 0.95)))),'k--'), 
plot(sort(longb),5 + y(indx_long),'k.','MarkerSize',10), hold
hold
set(gcf,'PaperUnits','centimeters')
set(gcf,'papersize',[15,10])
set(gcf,'paperposition',[0,0,15,10])
set(gcf,'units','centimeters')
set(gcf,'Position',[0,0,15,10])
xlabel('Longitude','interpreter','Latex','fontsize',17)
ylabel('Magnitude','interpreter','Latex','fontsize',17)
set(gca,'fontsize',13)
axs=axis
axis([minlong maxlong mu (mu + 3.0)])
hold off
box on
print(2, '-dpdf', 'expon_spatial_15_magnitude_vs_longitude.pdf')

% Sigma 
sigma0a = zeros(Lb,lx,m);
for ms = 1:m
    for j = 1:Lb
        [ms,j]
        z = (theta1a(j,2:Kgx,ms))';
        sigvec = TE*[theta1a(j,1,ms); (exp(theta1a(j,Kgx + 1,ms))*z)];
        sig = exp(X*sigvec);
        sigma0a(j,:,ms) = sig;  
    end
end

sigma0hata = zeros(1,lx);
sigma0mata = zeros(4*L,lx);
for k = 1:lx
sigma00a = sigma0a((burnin + 1):Lb,k,1);
sigma00b = sigma0a((burnin + 1):Lb,k,2);
sigma00c = sigma0a((burnin + 1):Lb,k,3);
sigma00d = sigma0a((burnin + 1):Lb,k,4);
sigma0mcmc = [sigma00a; sigma00b; sigma00c; sigma00d];
sigma0mata(:,k) = sigma0mcmc;  
sigma0hata(k) = median(sigma0mcmc);
end

% DIC and pd
Dhat3 = 0;
for i = 1:n
    Dhat3 = Dhat3 + 2*explike([sigma0hata(i)],y(i))
end
Dall3 = zeros(1,L*4);
for j = 1:(L*4)
    Dtemp = 0;
    for i = 1:n
        Dtemp = Dtemp + 2*explike([sigma0mata(j,i)],y(i));
    end   
    Dall3(j) = Dtemp;
end
Dave3 = mean(Dall3);
pd3 = Dave3 - Dhat3;
DIC3 = Dhat3 + 2*pd3

%%%%%%%%%%%%%% OLD STUFF  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(1), plot(long(indloc),lat(indloc),'r.')
% sigm1 = 0.7447;
% %y = -sigm1*log(1 - rand(2000,1));
% xi = -0.05;
% [gppar,gpparci] = gpfit(y)
% [exppar,expparci] = expfit(y)
% sigm = exppar; %mean(y);
% p = (1:n)/(n + 1);
% t = 0.01:0.01:2.5;
% ft = exp(-t/sigm);
% ft2 = (1 + gppar(1)*t/gppar(2)).^(-1/gppar(1));
% %plot(-log(1 - p),ysort,'k.')
% figure(2), plot(-sigm*log(1 - p),ysort,'*')
% hold
% plot(t,t,'k-')
% hold
% 
% figure(3), plot(t,1 - ft,'b-')
% hold
% plot(t,1 - ft2,'r-')
% plot(ysort,p,'k*')
% hold
% 
% % Y < sigma/(-xi)
% 
% % (-xi)*ymax = sigma
% 
% xiplot = -0.35:0.01:0;
% sigmaplot = (-xiplot)*max(y) ;
% 
% figure(4),plot(xiplot,sigmaplot,'b-')
% hold
% plot(gppar(1),gppar(2),'r*')
% plot(gpparci(:,1),gppar(2)*ones(2,1),'r-')
% plot(gppar(1)*ones(2,1),gpparci(:,2),'r-')
% hold
% 
% 
% plot(longb,sig,'k*')
% 
% [sort_long,indx_long] = sort(longb);
% 
% figure(5), plot(sort(y./(sig/(-xi))),'k*')
% figure(6), plot(sort(longb),5+(sig(indx_long)/(-xi)),'k-')
% hold
% plot(sort(longb),5+(sig(indx_long)*(2^(xi)-1)/(xi)),'r-')
% plot(sort(longb),5 + y(indx_long),'b*')
% hold
% xlabel('Longitude')
% ylabel('Magnitude')
% title('Magnitude data (blue), median (red), upper bound (black)')
% %print(6, '-dpdf', 'magnitude_vs_longitude.pdf')
% 
