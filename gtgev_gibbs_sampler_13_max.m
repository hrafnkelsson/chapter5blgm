% Bayesian modeling of spatially correlated extreme
% values using the multivariate Gaussian tail model.

% The parameters mu and tau are modeled with x, y and
% z = distance from coast. In addition tensor products of
% B-splines are used to model mu and tau in two dimensions. 
% The coefficients of the tensor products are given
% a Markov random field prior.

clear all
clear global

% The Gibbs sampler

% Load data on maximum temperature
load /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/data_matrix_max0

% Load data on distances between sites and the location of sites
load /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/dist02
load /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/xyloc02

% Load data giving the outline of Iceland
load /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/iceland.dat
%load /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/line_dist_sea.txt
sx = xyloc02(:,1);
sy = xyloc02(:,2);
alti = xyloc02(:,3);
dist_osea = xyloc02(:,4);

% Determine whether maxima (vswitch = 1) or minima (vswitch = -1) 
% is being modeled. 
vswitch = 1;
y = vswitch*data_matrix_max0;

% Clear a few variables
dist = dist02;
clear dist02 data_matrix_max0 xyloc02

% Determine the sign of gamma
sign_gamma = -1;

% Determine whether psi is a constant or a vector.
psi_constant = 1;

% Allocate the non-missing data
nonmiss1 = (y~=vswitch*(-9999));
[T,N] = size(y);
t0 = ceil(T/2)

% Create index matrices based on nonmiss1
C_cell = cell(N,1);
for i = 1:N
    C_cell{i}=find(nonmiss1(:,i));
end
i_t = zeros(T,N);
for i = 1:N
    for t = 1:T
        if nonmiss1(t,i) == 1
            i_t(t,i) = sum(nonmiss1(t,1:i));
        elseif nonmiss1(t,i) == 0
            i_t(t,i) = -9999;
        end
    end
end
A_cell = cell(T,1);
for t = 1:T
    A_cell{t} = find(nonmiss1(t,:));
end

% Number of observations per year.
nt = sum(nonmiss1');

% Number of observation per station.
ni = sum(nonmiss1);

% Total number of observations.
Ntot = sum(nt);

% Define the global variables
global T N Ntot t0 y nonmiss1 dist R_y
global R_y_inv_cell sign_gamma nt psi0xxx
global C_cell i_t A_cell X_mu X_phi XmmXm p_x_mu
global Z1 Z2 Z1mZ1 Z2mZ2 Nx1 Ny1 Nx2 Ny2
global X_tau XtmXt p_x_tau

global mu0 mu1 tau0 tau1
global Kx Ky Kmd eigC1 M1inv H1 eigC2 M2inv H2  
global alpha_phi_1_mu beta_phi_1_mu alpha_phi_2_mu beta_phi_2_mu 
global alpha_phi_1_tau beta_phi_1_tau alpha_phi_2_tau beta_phi_2_tau 
global eta_1_mu eta_1_tau 

global std_prop_mu std_prop_tau std_prop_psi
global mean_phi_y var_phi_y std_prop_phi_y
global mean_nu_y var_nu_y std_prop_nu_y

global mu_mu_0 s2_mu_0
global mu_tau_0 s2_tau_0
global mu_psi_0 s2_psi_0 std_prop_psi0
global mu_delta_0 s2_delta_0 std_prop_delta0
global p_psi s2_beta_psi mu_beta_psi
global s2_s2_psi mu_s2_psi
global s2_beta mu_beta s2_eta mu_eta
global s2_s2_mu nu_s2_mu s2_s2_tau nu_s2_tau
global s2_s2_mu_1 nu_s2_mu_1 s2_s2_mu_2 nu_s2_mu_2
global s2_kappa2_1_mu nu_kappa2_1_mu
global s2_kappa2_2_mu nu_kappa2_2_mu 
global s2_s2_tau_1 nu_s2_tau_1 s2_s2_tau_2 nu_s2_tau_2
global s2_kappa2_1_tau nu_kappa2_1_tau
global s2_kappa2_2_tau nu_kappa2_2_tau

% Compute an estimate of mu and tau using the Gumbel distribution.
% Fit parameters to individuals sites.
tau0 = zeros(N,1);
mu0 = zeros(N,1);
std_mu = zeros(N,1);
std_tau = zeros(N,1);
for i = 1:N
    index0 = find(nonmiss1(:,i));
    z0 = y(index0,i);
    [theta,theta_ci] = gevfit(z0);
    mu0(i) = theta(3);
    tau0(i) = log(theta(2));
    kappa0(i) = theta(1);
    std_mu(i) = (theta_ci(2,3)-theta_ci(1,3))/4;
    std_tau(i) = (theta_ci(2,2)-theta_ci(1,2))/4/theta(2);
    std_kappa(i) = (theta_ci(2,1)-theta_ci(1,1))/4;
end
index1 = find(1-isnan(std_mu));
index2 = find(isnan(std_mu));
nix2 = length(index2);
std_mu_average = mean(std_mu(index1));
std_mu(index2) = std_mu_average*ones(nix2,1);
std_tau_average = mean(std_tau(index1));
std_tau(index2) = std_tau_average*ones(nix2,1);
std_kappa_average = mean(std_kappa(index1));
std_kappa(index2) = std_kappa_average*ones(nix2,1);

tau0_all = 0.9;
tau0xxx = tau0_all*ones(N,1);
mu0xxx = zeros(N,1);
psi0xxx = 0.13;
global y00 tau0_all
for i = 1:N
    index0 = find(nonmiss1(:,i));
    y00 = y(index0,i);
    [theta,theta_ci] = evfit(-y00);
    mu0x(i) = -theta(1);
    tau0x(i) = log(theta(2));
    mu_tau = fminsearch('gev_gam_fixed_2',[mu0x(i) tau0x(i)]);
    mu0xxx(i) = mu_tau(1);
    tau0xxx(i) = mu_tau(2);
end

u00 = [];
for i = 1:N
    index0 = find(nonmiss1(:,i));
    y00 = y(index0,i);
    u00 = [u00 (y00'-mu0xxx(i))/exp(tau0xxx(i))];
end
[theta_u,theta_ci_u] = gevfit(u00);
std_psi0 = (theta_ci_u(2,1)-theta_ci_u(1,1))/4;

% MCMC: Metropolis--Hastings steps within the Gibbs sampler

% Define the covariance matrices for each year
R_y_inv_cell = cell(T,1);
for j = 1:T
    Ry_inv_cell{j} = zeros(nt(j),nt(j));
end

% Define the fixed effects regression matrices for mu and tau
load /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/Lon.out
load /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/Lat.out
med_lon = median(Lon);
min_lon = min(Lon);
max_lon = max(Lon);
med_lat = median(Lat);
min_lat = min(Lat);
max_lat = max(Lat);
clear Lon Lat
% Skip longitude
X_mu=zeros(N,3);
X_mu(:,1) = ones(N,1);
%X_mu(:,2) = sy - med_lat; 
X_mu(:,2) = alti/100;
X_mu(:,3) = dist_osea/100;

p_x_mu = size(X_mu,2);
XmmXm = X_mu'*X_mu;

% Create a specific fixed effect matirx for tau. Call it X_tau
% Skip longitude and latitude
X_tau = zeros(N,1);
X_tau(:,1) = ones(N,1);
%X_tau(:,2) = sx - med_lon; 
%X_tau(:,3) = sy - med_lat; 

p_x_tau = size(X_tau,2);
XtmXt = X_tau'*X_tau;

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

% The number of equally spaced interior knots.
%kx = 56;
%ky = 38;

kx = 56;
ky = 38;

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
x_bsp = (sx - min_lon)/(max_lon - min_lon);
y_bsp = (sy - min_lat)/(max_lat - min_lat);
lx = length(x_bsp);
ly = length(y_bsp);

% Compute the x-splines and the y-splines.
X = spline_functions(x_bsp,tau_x,dx,kx,M);
Y = spline_functions(y_bsp,tau_y,dy,ky,M);

% Create the matrix Z1 with the tensor product of X and Y.
Z1 = zeros(N,Nx1*Ny1);
for i = 1:Ny1
    for j = 1:Nx1
        Z1(:,Nx1*(i-1)+j) = X(:,j).*Y(:,i);
    end
end

%-------------------------------------------------------------------------

% Z2 : a matrix with high-dimensional B-spline tensor-products

% The number of equally spaced interior knots.
kx = 26;
ky = 17;

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
x_bsp = (sx - min_lon)/(max_lon - min_lon);
y_bsp = (sy - min_lat)/(max_lat - min_lat);
lx = length(x_bsp);
ly = length(y_bsp);

% Compute the x-splines and the y-splines.
X = spline_functions(x_bsp,tau_x,dx,kx,M);
Y = spline_functions(y_bsp,tau_y,dy,ky,M);

% Create the matrix Z2 with the tensor product of X and Y.
Z2 = zeros(N,Nx2*Ny2);
for i = 1:Ny2
    for j = 1:Nx2
        Z2(:,Nx2*(i-1)+j) = X(:,j).*Y(:,i);
    end
end

%-------------------------------------------------------------------------

Z1mZ1 = Z1'*Z1;
Z2mZ2 = Z2'*Z2;

% Define the regression matrices for phi
ace = ones(N,1);
if psi_constant ~= 1
    X_phi = ones(N,1);
    p_phi = size(X_phi,2);
end

% Define the matrices H1, H2, C1, C2, M1 and M2.

%-------------------------------------------------------------------------

% H1, C1 and M1
Kgx = Nx1;
Kgy = Ny1;
Rx = zeros(1,Kgx);
Rx(2) = 1;
Rx = toeplitz(Rx);
Ry = zeros(1,Kgy);
Ry(2) = 1;
Ry = toeplitz(Ry);
H1 = kron(eye(Kgy),Rx)+kron(Ry,eye(Kgx));
H1 = sparse(H1);
M1 = sum(H1');
M1 = sparse(diag(1./M1));

for i = 1:Kgx*Kgy    
    C1(i,:) = H1(i,:)*M1(i,i);
    M1inv(i,i) = 1/M1(i,i);
end
I1 = sparse(eye(Kgx*Kgy));
M1inv = sparse(M1inv);
'eigenvalues1'
% eigC1 = real(eig(full(C1)));
% save /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/eigC1_60_42 eigC1
load /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/eigC1_60_42

%-------------------------------------------------------------------------

% H2, C2 and M2
Kgx = Nx2;
Kgy = Ny2;
Rx = zeros(1,Kgx);
Rx(2) = 1;
Rx = toeplitz(Rx);
Ry = zeros(1,Kgy);
Ry(2) = 1;
Ry = toeplitz(Ry);
H2 = kron(eye(Kgy),Rx)+kron(Ry,eye(Kgx));
H2 = sparse(H2);
M2 = sum(H2');
M2 = sparse(diag(1./M2));

for i = 1:Kgx*Kgy
    C2(i,:) = H2(i,:)*M2(i,i);
    M2inv(i,i) = 1/M2(i,i);
end
I2 = sparse(eye(Kgx*Kgy));
M2inv = sparse(M2inv);
'eigenvalues2' 
% eigC2 = real(eig(full(C2)));
% save /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/eigC2_30_21 eigC2
load /home/birgirhr/herdubreid/hydrology_meteorology/meteorology_gogn/eigC2_30_21

%-------------------------------------------------------------------------

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

