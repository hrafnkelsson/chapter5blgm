% Bayesian modeling of spatially varying earthquake magnitudes
% that are model with the generalized Pareto distribution (GPD).
% Assume sigma of the GPD varies spatially while xi is fixed
% and the same for the whole area.

% Tensor products of B-splines are used to model theta=log(sigma)
% in two dimensions.
% The coefficients of the tensor products are given
% a Markov random field prior.

clear all
clear global

global Q1 Z1 K y lam1 lam2

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
xi = -0.16;
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

% The spatial dimensions of the rectangular area.
min_lon = -22.5;
max_lon = -19.8;
min_lat = 63.75;
max_lat = 64.15;
% lon 132 km
% lat 33 km
% % Define the random regression matrices Z1, Z2 and Z3 for mu and tau
% Kx = 20;    % The number of knots for the x-direction main effect.
% Ky = 20;    % The number of knots for the y-direction main effect.
% Kmd = 12;   % The number of knots for the min-direction main effect.
% knots_x = ([1:Kx])/(Kx+1);
% knots_x = prctile(x_scaled,knots_x*100);
% knots_y = ([1:Ky])/(Ky+1);
% knots_y = prctile(y_scaled,knots_y*100);
% knots_md = ([1:Kmd])/(Kmd+1);
% knots_md = prctile(md_scaled,knots_md*100);

%-------------------------------------------------------------------------

% Z1 : a matrix with low dimensional B-spline tensor products
% Z1 is used to infer the unknown B-spline coefficients.

% The number of equally spaced interior knots.
kx = 24 - 4;
ky = 6 - 4;
%kx = 16;
%ky = 6;


% Delta x and Delta y.
dx = 1/(kx+1);
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
sy_temp = min_lat + (max_lat - min_lat)*(0:(1/25):1);
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

%-------------------------------------------------------------------------

% Define the matrix Q1.

%-------------------------------------------------------------------------

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
Qchol = chol(Q1);
% diffQ=Q1 - Qchol'*Qchol;
% sum(sum(diffQ))

% Simple inference - the coefficients:
lam1 = 1/0.5;
lam2 = 1/(0.5/Nx1);
v1 = 0.22 - 0.02*((0:(Nx1 - 1))/(Nx1 - 1));
tst = kron(ones(Ny1,1),v1');
startvec = [tst; 0.47; log(0.03)];
postdens_dim2(startvec)
[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim2(theta),startvec);
postdens_dim2(t_m)

% [t_m_2,~,~,~,~,H] = fminunc(@(theta) postdens_dim2(theta),t_m);
% postdens_dim2(t_m_2)

stop

% sigma_nu = 0.5;
% Qpost = ((1/sigma_nu^2)*Q1 + Z1'*Qdata0*Z1);   %
% mupost = Qpost\(Z1'*Qdata0*muvec0);             %
mupostmatrix = reshape(t_m(1:K),Nx1,Ny1);
xcont = min_lon + (max_lon - min_lon)*(0:1:(Nx1 - 1))/(Nx1 - 1);
ycont = min_lat + (max_lat - min_lat)*(0:1:(Ny1 - 1))/(Ny1 - 1);
figure(1)
surf(xcont,ycont,(mupostmatrix)')
view([0 90])
colorbar

% Simple inference - eta = Z*nu:
etapost = Z2*t_m(1:K);
etapostmatrix = reshape(etapost,lx_temp,ly_temp);
xcont = min_lon + (max_lon - min_lon)*(0:1:(lx_temp - 1))/(lx_temp - 1);
ycont = min_lat + (max_lat - min_lat)*(0:1:(ly_temp - 1))/(ly_temp - 1);
figure(2)
surf(xcont,ycont,(etapostmatrix)'); 
view([0 90])
colorbar
axis([min_lon max_lon min_lat max_lat])
xlabel('Longitude')
ylabel('Latitude')
title('Sigma')
print(2, '-dpdf', 'sigma_oef_dim2.pdf')


% The true field
%mufcn2 = 0.5 + 0.15*(x_bsp_temp - 1) + 0.15*(y_bsp_temp - 0.5) + 0.15*(x_bsp_temp - 1).^2 - 0.10*(y_bsp_temp - 0.5).^2;
mufcn2 = 0.5 + 0.025*cos(8*x_bsp_temp - 0.3) + 0.025*cos(8*y_bsp_temp - 0.2);
mufcn2matrix = reshape(mufcn2,lx_temp,ly_temp);
figure(3)
surf(sx_temp,sy_temp,(mufcn2matrix)')
view([0 90])
colorbar

% The difference between estimate and true 
figure(4)
surf(sx_temp,sy_temp,(lowb + diffb*(1+exp(-etapostmatrix)).^(-1))'-(mufcn2matrix)')
view([0 90])
colorbar

ratio_median = median(median((lowb + diffb*(1+exp(-etapostmatrix)).^(-1))'./(mufcn2matrix)'))
log_diff_median = median(median(log(lowb + diffb*(1+exp(-etapostmatrix)).^(-1))'-log(mufcn2matrix)'))

% Estimate the posterior value for sigma_nu





%%%%%%%%%%%%%%%%%%%%%
STOP
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% old stuff for MCMC ...

% Define the parameters in the hyper-priors
% phi_y and nu_y
mean_phi_y = -7;
var_phi_y = 4^2;
mean_nu_y = -1;
var_nu_y = 1.0^2;

% mu_mu_0, s2_mu_0
mu_mu_0 = 0;
s2_mu_0 = 1000^2;

% mu_tau_0, s2_tau_0
mu_tau_0 = 0;
s2_tau_0 = 1000^2;

% mu_psi_0, s2_psi_0
mu_psi_0 = -1.6;
s2_psi_0 = 2^2;

% mu_delta_0, s2_delta_0
mu_delta_0 = 0;
s2_delta_0 = 100^2;

%s2_beta_mu mu_beta_mu, s2_eta_mu mu_eta_mu, 
s2_beta = 1000000;
mu_beta = 0;
s2_eta = 1000000;
mu_eta = 0;

%s2_s2_mu nu_s2_mu
s2_s2_mu = 1;
nu_s2_mu = 10^(-12);

%s2_s2_mu_j nu_s2_mu_j, j=1,2,3
s2_s2_mu_1 = 1;
nu_s2_mu_1 = 10^(-12);
s2_kappa2_1_mu = 1;
nu_kappa2_1_mu = 10^(-12);

%s2_s2_tau nu_s2_tau
s2_s2_tau = 1;
nu_s2_tau = 10^(-12);

%s2_s2_tau_j nu_s2_tau_j, j=1,2,3
s2_s2_tau_1 = 1;
nu_s2_tau_1 = 10^(-12);
s2_kappa2_1_tau = 1;
nu_kappa2_1_tau = 10^(-12);

%alpha_phi beta_phi
alpha_phi_1_mu = 100;
beta_phi_1_mu = 0.5;

alpha_phi_1_tau = 100;
beta_phi_1_tau = 0.5;

if psi_constant ~= 1
    % beta_psi, s2_psi, phi_psi and nu_psi
    s2_beta_psi = 1000^2;
    mu_beta_psi = 0;
    s2_s2_psi = 0.01;
    mu_s2_psi = 0.01;
end

% The size of the MCMC chain
L = 8000;
burnin = 2000;
runtime = burnin+L;
mc = 4;

% L = 4000;
% burnin = 1500;
% runtime = burnin+L;
% mc = 3;

% Define the variables
mu_gen = zeros(N,runtime,mc);
delta0_gen = zeros(runtime,mc);
tau_gen = zeros(N,runtime,mc);
if psi_constant == 1
    psi0_gen = zeros(runtime,mc);
else
    psi_gen = zeros(N,runtime,mc);
end
phi_y_gen = zeros(runtime,mc);
nu_y_gen = zeros(runtime,mc);
if psi_constant ~= 1
    beta_psi_gen = zeros(runtime,mc);
    s2_psi_gen = zeros(runtime,mc);
end
D = zeros(runtime,mc);
beta_gen = zeros(p_x_mu,runtime,mc);
eta_gen = zeros(p_x_tau,runtime,mc);
s2_mu_gen = zeros(runtime,mc);
s2_tau_gen = zeros(runtime,mc);
a1_gen = zeros(Nx1*Ny1,runtime,mc);
b1_gen = zeros(Nx2*Ny2,runtime,mc);
kappa2_1_mu_gen = zeros(runtime,mc);
phi_1_mu_gen = zeros(runtime,mc);
kappa2_1_tau_gen = zeros(runtime,mc);
phi_1_tau_gen = zeros(runtime,mc);

% Random initial values
if psi_constant == 1
    psi0_gen(1,1:mc) = log(psi0xxx)*ones(1,mc)+0*randn(1,mc)./sqrt(N*T);
elseif psi_constant ~= 1
    psi_gen(:,1,1:mc) = log(kron(psi0xxx*ones(N,1),ones(1,mc)).*gamrnd(N*T/2,2/T/N,N,mc));
end
tau_gen(:,1,1:mc) = log(kron(exp(tau0xxx),ones(1,mc)).*gamrnd(N*T/2,2/T/N,N,mc));
mu_gen(:,1,1:mc) = kron(mu0xxx,ones(1,mc));
for m = 1:mc
    for i = 1:N
        i
        fyi = 0;
        while fyi == 0
            v = mu_gen(i,1,m)+0.01*randn(1,1)*std_mu(i);
            if psi_constant == 1
                fyi = min(gevpdf(y(find(nonmiss1(:,i)==1),i),sign_gamma*exp(psi0_gen(1,m)),exp(tau_gen(i,1,m)),v));
            elseif psi_constant ~= 1
                fyi = min(gevpdf(y(find(nonmiss1(:,i)==1),i),sign_gamma*exp(psi_gen(i,1,m)),exp(tau_gen(i,1,m)),v));
            end
        end
        mu_gen(i,1,m) = v;
    end
end
psi0_gen = kron(log(0.14),ones(1,mc)).*gamrnd(N*T/2,2/T/N,1,mc);
delta0_gen(1,1:mc) = 0.025 + 0.005*randn(1,mc);    
phi_y_gen(1,1:mc) = kron(log(0.0010),ones(1,mc)).*gamrnd(N*T/2,2/T/N,1,mc);
nu_y_gen(1,1:mc) = kron(log(0.15),ones(1,mc)).*ones(1,mc); %  .*gamrnd(N*T/2,2/T/N,1,mc);
if psi_constant~=1
    beta_psi_gen(1,1:mc)=kron((-0.7)*ones(p_psi,1),ones(1,mc))+randn(p_psi,mc).*kron(mean(exp(tau0)),ones(1,mc))/T;
    s2_psi_gen(1,1:mc)=kron(0.2,ones(1,mc)).*gamrnd(N*T/2,2/T/N,1,mc);
end
s2_mu_gen(1,1:mc) = kron(1,0.8*ones(1,mc)).*gamrnd(N*T/2,2/T/N,1,mc);
s2_tau_gen(1,1:mc) = kron(0.004,ones(1,mc)).*gamrnd(N*T/2,2/T/N,1,mc);
kappa2_1_mu_gen(1,:) = 0.5*ones(1,mc)*10^(1);
phi_1_mu_gen(1,:) = ones(1,mc)*0.985;
kappa2_1_tau_gen(1,:) = 3*ones(1,mc)*10^(-2);
phi_1_tau_gen(1,:) = ones(1,mc)*0.985;

for k = 1:mc
    beta_gen(:,1,k) = inv(X_mu'*X_mu)*X_mu'*mu0xxx + 0.01*randn(p_x_mu,1);    
    a1_gen(:,1,k) = sqrt(kappa2_1_mu_gen(1,k))*randn(Nx1*Ny1,1);
    eta_gen(:,1,k) = inv(X_tau'*X_tau)*X_tau'*tau0xxx + 0.01*randn(p_x_tau,1);    
    b1_gen(:,1,k) = sqrt(kappa2_1_tau_gen(1,k))*randn(Nx2*Ny2,1);
end

% Proposal variances for mu, tau and psi.
std_prop_mu = std_mu*2.2;
std_prop_tau = std_tau*2.2;
eta_1_mu = 0.1;
eta_1_tau = 0.08;

std_prop_phi_y = 0.45;
std_prop_nu_y = 0.07;
std_prop_delta0 = 0.03;
if psi_constant == 1
    std_prop_psi0 = std_psi0*14.0;
elseif psi_constant ~= 1
    std_prop_psi = 2*mean(std_kappa.*sqrt(ni))./ni';
end

if psi_constant == 1

    % Start the Gibbs sampler which uses Metropolis-Hastings steps.
    for m = 1:mc                                                        
        tic
        % Initialize the matrix Ry
        R_y = matcorr0([exp(phi_y_gen(1,m)) exp(nu_y_gen(1,m))]);
        for t = 1:T
            R_y_inv_cell{t} = inv(R_y(A_cell{t},A_cell{t}));
        end
        mu0 = X_mu*beta_gen(:,1,m); mu1 = Z1*a1_gen(:,1,m); 
        tau0 = X_tau*eta_gen(:,1,m); tau1 = Z2*b1_gen(:,1,m); 
        for j = 2:runtime
            [m j]
            %[m j 1]
            %mu_gen(:,j,m)=mu0xxx;            
            for i = 1:N
                mu_gen(i,j,m) = rgen_mu([mu_gen(1:i-1,j,m)' mu_gen(i:N,j-1,m)']',tau_gen(:,j-1,m), ...
                    psi0_gen(j-1,m),delta0_gen(j-1,m),s2_mu_gen(j-1,m),i);
            end
            %[m j 2]
            %tau_gen(:,j,m)=tau0xxx;
            for i = 1:N
                tau_gen(i,j,m)=rgen_tau(mu_gen(:,j,m),[tau_gen(1:i-1,j,m)' tau_gen(i:N,j-1,m)']', ...
                    psi0_gen(j-1,m),delta0_gen(j-1,m),s2_tau_gen(j-1,m),i);
            end   
            psi0_gen(j,m) = rgen_psi0(mu_gen(:,j,m),tau_gen(:,j,m),psi0_gen(j-1,m),delta0_gen(j-1,m));
            delta0_gen(j,m) = rgen_delta0(mu_gen(:,j,m),tau_gen(:,j,m),psi0_gen(j,m),delta0_gen(j-1,m));                           
            phi_y_gen(j,m) = rgen_phi_y(mu_gen(:,j,m),tau_gen(:,j,m),psi0_gen(j,m),delta0_gen(j,m),phi_y_gen(j-1,m),nu_y_gen(j-1,m));            
            nu_y_gen(j,m) = rgen_nu_y(mu_gen(:,j,m),tau_gen(:,j,m),psi0_gen(j,m),delta0_gen(j,m),phi_y_gen(j,m),nu_y_gen(j-1,m));            
            beta_gen(:,j,m) = rgen_beta(mu_gen(:,j,m),s2_mu_gen(j-1,m)); mu0 = X_mu*beta_gen(:,j,m);
            eta_gen(:,j,m) = rgen_eta(tau_gen(:,j,m),s2_tau_gen(j-1,m));  tau0 = X_tau*eta_gen(:,j,m);
            s2_mu_gen(j,m) = rgen_s2_mu(mu_gen(:,j,m));
            % s2_mu_gen(j,m) = 0.01;
            s2_tau_gen(j,m) = rgen_s2_tau(tau_gen(:,j,m));
            % s2_tau_gen(j,m) = 0.0001;
            a1_gen(:,j,m) = rgen_a1(mu_gen(:,j,m),s2_mu_gen(j,m),kappa2_1_mu_gen(j-1,m),phi_1_mu_gen(j-1,m)); mu1 = Z1*a1_gen(:,j,m);
            kappa2_1_mu_gen(j,m) = rgen_kappa2_1_mu(a1_gen(:,j,m),phi_1_mu_gen(j-1,m));
            phi_1_mu_gen(j,m) = rgen_phi_1_mu(a1_gen(:,j,m),kappa2_1_mu_gen(j,m),phi_1_mu_gen(j-1,m));
            b1_gen(:,j,m) = rgen_b1(tau_gen(:,j,m),s2_tau_gen(j,m),kappa2_1_tau_gen(j-1,m),phi_1_tau_gen(j-1,m)); tau1 = Z2*b1_gen(:,j,m);
            kappa2_1_tau_gen(j,m) = rgen_kappa2_1_tau(b1_gen(:,j,m),phi_1_tau_gen(j-1,m));
            phi_1_tau_gen(j,m) = rgen_phi_1_tau(b1_gen(:,j,m),kappa2_1_tau_gen(j,m),phi_1_tau_gen(j-1,m));
        end
    end
    toc

elseif psi_constant ~= 1

    % Start the Gibbs sampler which uses Metropolis-Hastings steps.
    for m = 1:mc
        % Initialize the matrix Ry
        R_y = matcorr0([exp(phi_y_gen(1,m)) exp(nu_y_gen(1,m))]);
        %R_mu_inv = inv(matcorr0([exp(phi_mu_gen(1,m)) exp(nu_mu_gen(1,m))]));
        %R_psi_inv=inv(matcorr0([phi_psi_gen(1,m) nu_psi_gen(1,m)]));
        for t = 1:T
            R_y_inv_cell{t} = inv(R_y(A_cell{t},A_cell{t}));
        end
        for j = 2:runtime
            [m j]
            [m j 1]
            %mu_gen(:,j,m) = mu0xxx;
            for i = 1:N
                mu_gen(i,j,m) = rgen_mu([mu_gen(1:i-1,j,m)' mu_gen(i:N,j-1,m)']',tau_gen(:,j-1,m),...
                    psi_gen(:,j-1,m),phi_y_gen(j-1,m),nu_y_gen(j-1,m),...
                    beta_mu_gen(j-1,m),s2_mu_gen(j-1,m),phi_mu_gen(j-1,m),nu_mu_gen(j-1,m),i);
            end
            [m j 2]
            %tau_gen(:,j,m)=tau0xxx;
            for i = 1:N
                tau_gen(i,j,m)=rgen_tau(mu_gen(:,j,m),[tau_gen(1:i-1,j,m)' tau_gen(i:N,j-1,m)']', ...
                    psi_gen(:,j-1,m),phi_y_gen(j-1,m),nu_y_gen(j-1,m),...
                    beta_tau_gen(j-1,m),s2_tau_gen(j-1,m),i);
            end
            [m j 3]%i
            for i = 1:N
                %[m j 3 i]%i
                psi_gen(i,j,m) = rgen_psi(mu_gen(:,j,m),tau_gen(:,j,m), ...
                    [psi_gen(1:i-1,j,m)' psi_gen(i:N,j-1,m)']',phi_y_gen(j-1,m),nu_y_gen(j-1,m),...
                    beta_psi_gen(j-1,m),s2_psi_gen(j-1,m),i);
            end
            phi_y_gen(j,m) = rgen_phi_y(mu_gen(:,j,m),tau_gen(:,j,m),psi_gen(:,j,m),phi_y_gen(j-1,m),nu_y_gen(j-1,m));
            %phi_y_gen(j,m) = log(1/10);
            nu_y_gen(j,m) = rgen_nu_y(mu_gen(:,j,m),tau_gen(:,j,m),psi_gen(:,j,m),phi_y_gen(j,m),nu_y_gen(j-1,m));
            %nu_y_gen(j,m) = log(0.5);
            beta_psi_gen(j,m) = rgen_beta_psi_2(psi_gen(:,j,m),s2_psi_gen(j-1,m));
            s2_psi_gen(j,m) = rgen_s2_psi_2(psi_gen(:,j,m),beta_psi_gen(j,m));
            D(j,m) = D_y_theta(mu_gen(:,j,m),tau_gen(:,j,m),psi_gen(j,m),phi_y_gen(j,m),nu_y_gen(j,m));
        end
    end

end  % end if statement about psi

% % DIC computations
% for m = 1:mc
%     for j = 1:runtime
%         [m j]
%         R_y = matcorr0([exp(phi_y_gen(j,m)) exp(nu_y_gen(j,m))]);
%         for t = 1:T
%             R_y_inv_cell{t} = inv(R_y(A_cell{t},A_cell{t}));
%         end
%         D(j,m) = D_y_theta(mu_gen(:,j,m),tau_gen(:,j,m),psi0_gen(j,m),phi_y_gen(j,m),nu_y_gen(j,m));
%     end
% end

save gtgev_mu_60_42_tau_30_21_10000_max
% 
% gtgev_9_prediction

stop1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if psi_constant == 1
    mc = 1
    mu_median = zeros(N,1);
    tau_median = zeros(N,1);
    for k = 1:N
        mu_median(k) = median(reshape(mu_gen(k,:,1:mc),[L*mc 1]));
        tau_median(k) = median(reshape(tau_gen(k,:,1:mc),[L*mc 1]));
    end
    psi0_median = median(reshape(psi0_gen(:,1:mc),[L*mc 1]));
    phi_y_median = median(reshape(phi_y_gen(:,1:mc),[L*mc 1]));
    nu_y_median = median(reshape(nu_y_gen(:,1:mc),[L*mc 1]));
    R_y = matcorr0([exp(phi_y_median) exp(nu_y_median)]);
    for t = 1:T
        R_y_inv_cell{t} = inv(R_y(A_cell{t},A_cell{t}));
    end
    D_hat_1 = D_y_theta(mu_median,tau_median,psi0_median,phi_y_median,nu_y_median)
    D_ave = mean(mean(D([11:L],1:mc)))
    DIC_1 = 2*D_ave-D_hat_1
    p_D_1 = D_ave-D_hat_1
elseif psi_constant ~= 1
    mu_median = zeros(N,1);
    tau_median = zeros(N,1);
    psi_median = zeros(N,1);
    for k = 1:N
        mu_median(k) = median(reshape(mu_gen(k,:,:),[L*mc 1]));
        tau_median(k) = median(reshape(tau_gen(k,:,:),[L*mc 1]));
        psi_median(k) = median(reshape(psi_gen(k,:,:),[L*mc 1]));
    end
    phi_y_median = median(reshape(phi_y_gen(:,:),[L*mc 1]));
    nu_y_median = median(reshape(nu_y_gen(:,:),[L*mc 1]));
    R_y = matcorr0([exp(phi_y_median) exp(nu_y_median)]);
    for t = 1:T
        R_y_inv_cell{t} = inv(R_y(A_cell{t},A_cell{t}));
    end
    D_hat_2 = D_y_theta(mu_median,tau_median,psi_median,phi_y_median,nu_y_median)
    D_ave = mean(mean(D([11:L],1:mc)))
    DIC_2 = 2*D_ave-D_hat_2
    p_D_2 = D_ave-D_hat_2

end

%ar_phi_y = sum(diff(phi_y_gen(:,1))~=0)/L
%ar_nu_y = sum(diff(nu_y_gen(:,1))~=0)/L
%ar_phi_mu = sum(diff(phi_mu_gen(:,1))~=0)/L
%ar_nu_mu = sum(diff(nu_mu_gen(:,1))~=0)/L
%ar_psi_0 = sum(diff(psi0_gen(:,1))~=0)/L

%save run_gtgev_model_1_b

mu_med = (median(mu_gen(:,501:1000,1)'))';
a1_med = (median(a1_gen(:,501:1000,1)'))';
beta_med = (median(beta_gen(:,501:1000,1)'))';
mu_pred = X_mu*beta_med + Z1*a1_med;
plot(sy,mu_med - mu_pred,'.')

