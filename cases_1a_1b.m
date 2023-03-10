% Analysis for Case 1a and Case 1b

% PSHA constant spatial model for the exponential distribution
% and the generalized Pareto distribution

clear all
clear global
close all

global y lam1
lam1 = 1/0.5;

% load the data 
[magn,long,lat] = readvars('ICEL_Mc4_SW.xlsx','Range','A2:C184');
plot(long,lat,'k.')
indloc = find((lat < 64.2).*(lat > 63.7).*(magn > 4.999));
plot(long(indloc),lat(indloc),'r.')

magnb = magn(indloc);
longb = long(indloc);
latb = lat(indloc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1a
y = magnb - 4.99;
ysort = sort(y);
n = length(y);
[exppar,expparci] = expfit(y);
ybar = mean(y);
a0 = 1;
b0 = lam1;
L = 4*25000;
beta = gamrnd((a0 + n),1/(b0 + n*ybar),L,1); 
betahat = median(beta);
Dhat = 2*explike(1/betahat,y);
Dall = zeros(1,L);
for j = 1:L
    Dall(j) = 2*explike(1/beta(j),y);
end
Dave = mean(Dall);
pd = Dave - Dhat;
DIC = Dhat + 2*pd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1b
[gppar,gpparci] = gpfit(y)
[negloglik,asymcov] = gplike(gppar,y)

theta0 = [log(0.5); 0.9];
pipost = postdens_dim0(theta0)

% Find the posterior mode
[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim0(theta),[log(0.48); 0.8]);
% non Z
%[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim1(theta),[log(0.48); randn(Kgx - 1,1); 0.5; log(0.001)]);
postdens_dim0(t_m)

% Metropolis step for the unknown parameters
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
    logpoststar = postdens_dim0(thetastar');
    if logpoststar == Inf 
         theta1(j,:,ms) = theta1(j - 1,:,ms);
    else
        logr = - logpoststar + postdens_dim0((theta1(j - 1,:,ms))');    % Step 2
        if log(rand(1)) < logr
            theta1(j,:,ms) = thetastar;  
        else 
            theta1(j,:,ms) = theta1(j - 1,:,ms);
        end % Step 3.
    end
end
end

pnr = 1;
plot(theta1(:,pnr,1),'b-'), hold,
plot(theta1(:,pnr,2),'k-')
plot(theta1(:,pnr,3),'r-')
plot(theta1(:,pnr,4),'g-'), hold
[gr,neff,ac,ar] = gpar(reshape(theta1(burnin+1:Lb,pnr,:),L,m));
gr

figure(1)
plot(exp(theta1(:,pnr,1)),'b-'), hold,
plot(exp(theta1(:,pnr,2)),'k-')
plot(exp(theta1(:,pnr,3)),'r-')
plot(exp(theta1(:,pnr,4)),'g-'), hold

pnr = 2;
xiref = 0.75;
figure(2)
% %xi = (-xiref + xiref*exp(theta1(:,pnr,1))/(1 + exp(theta1(:,pnr,1))));
plot((-xiref + xiref*exp(theta1(:,pnr,1))./(1 + exp(theta1(:,pnr,1)))),'b-'), hold,
plot((-xiref + xiref*exp(theta1(:,pnr,2))./(1 + exp(theta1(:,pnr,2)))),'k-')
plot((-xiref + xiref*exp(theta1(:,pnr,3))./(1 + exp(theta1(:,pnr,3)))),'r-')
plot((-0 + xiref*exp(theta1(:,pnr,4))./(1 + exp(theta1(:,pnr,4)))),'g-'), hold
psi2 = [theta1(:,pnr,1); theta1(:,pnr,2); theta1(:,pnr,3); theta1(:,pnr,4)];
figure(3), hist(psi2,100)

pnr = 2;
xi2a = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,1))./(1 + exp(theta1((burnin + 1):Lb,pnr,1)));
xi2b = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,2))./(1 + exp(theta1((burnin + 1):Lb,pnr,2)));
xi2c = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,3))./(1 + exp(theta1((burnin + 1):Lb,pnr,3)));
xi2d = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,4))./(1 + exp(theta1((burnin + 1):Lb,pnr,4)));
xi2mcmc = [xi2a; xi2b; xi2c; xi2d];
figure(4), hist(xi2mcmc,100)
xihat = median(xi2mcmc)

pnr = 1;
sigma2a = exp(theta1((burnin + 1):Lb,pnr,1));
sigma2b = exp(theta1((burnin + 1):Lb,pnr,2));
sigma2c = exp(theta1((burnin + 1):Lb,pnr,3));
sigma2d = exp(theta1((burnin + 1):Lb,pnr,4));
sigma2mcmc = [sigma2a; sigma2b; sigma2c; sigma2d];
figure(5), hist(sigma2mcmc,100)
sigmahat = median(sigma2mcmc)

% DIC and pd
Dhat2 = 2*gplike([xihat, sigmahat],y);
Dall2 = zeros(1,L*4);
for j = 1:(L*4)
    Dall2(j) = 2*gplike([xi2mcmc(j), sigma2mcmc(j)],y);
end
Dave2 = mean(Dall2);
pd2 = Dave2 - Dhat2;
DIC2 = Dhat2 + 2*pd2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figures
% Generalized Pareto distribution, constant sigma
minlong = -23;
maxlong = -19;
mu = 5;
[sort_long,indx_long] = sort(longb);

figure(1), plot([minlong maxlong],ones(1,2)*(mu+(sigmahat/(-xihat))),'k-')
hold
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.05)^(-xihat) - 1)/(xihat))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.25)^(-xihat) - 1)/(xihat))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.5)^(-xihat) - 1)/(xihat))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.75)^(-xihat) - 1)/(xihat))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.95)^(-xihat) - 1)/(xihat))),'k--')
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
print(1, '-dpdf', 'gp_const_magnitude_vs_longitude.pdf')

% Exponential distribution, constant sigma
minlong = -23;
maxlong = -19;
mu = 5;
[sort_long,indx_long] = sort(longb);

figure(2), hold
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.05)))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.25)))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.5)))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.75)))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.95)))),'k--')
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
print(2, '-dpdf', 'expon_const_magnitude_vs_longitude.pdf')

stop


%%%%%% OLD STUFF   %%%%%%%%%%%%%%%%%%%%%%%


% sigm1 = 0.7447;
% %y = -sigm1*log(1 - rand(2000,1));
% xi = -0.05;
% y = magnb - 4.99;
% ysort = sort(y);
% n = length(y);
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
