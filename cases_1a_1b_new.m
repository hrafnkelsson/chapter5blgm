% Analysis for Case 1a and Case 1b

% PSHA constant spatial model for the exponential distribution
% and the generalized Pareto distribution

clear all
clear global
close all

global y lam1 xiref alp bet 
lam1 = 1/0.5;
alp = 1;
bet = 1;
xiref = 0.75;

% load the data 
[magn,long,lat] = readvars('ICEL_Mc4_SW.xlsx','Range','A2:C184');
plot(long,lat,'k.')
indloc = find((lat < 64.025).*(lat > 63.585).*(magn > 4.499));

% long-min: -22.90.  long-max: -19.00.   
% lat-min:   63.59.  lat-max:   64.95.     

plot(long(indloc),lat(indloc),'r.')

magnb = magn(indloc);
longb = long(indloc);
latb = lat(indloc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1a
y = magnb - 4.49;
ysort = sort(y);
n = length(y);
[exppar,expparci] = expfit(y);
ybar = mean(y);
a0 = 1;
b0 = lam1;
L = 4*100000;
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

theta0 = [log(0.8); 0.9];
pipost = postdens_dim0(theta0)

% Find the posterior mode
[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim0(theta),[log(0.8); 0.9]);
% non Z
%[t_m,~,~,~,~,H] = fminunc(@(theta) postdens_dim1(theta),[log(0.48); randn(Kgx - 1,1); 0.5; log(0.001)]);
postdens_dim0(t_m)

% Metropolis step for the unknown parameters
K = length(t_m);
uc = 2.38*sqrt(1/K);
L = 100000;
burnin = 5000;
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
% plot(theta1(:,pnr,1),'b-'), hold,
% plot(theta1(:,pnr,2),'k-')
% plot(theta1(:,pnr,3),'r-')
% plot(theta1(:,pnr,4),'g-'), hold
[gr,neff,ac,ar] = gpar(reshape(theta1(burnin+1:Lb,pnr,:),L,m));
gr

% figure(1)
% plot(exp(theta1(:,pnr,1)),'b-'), hold,
% plot(exp(theta1(:,pnr,2)),'k-')
% plot(exp(theta1(:,pnr,3)),'r-')
% plot(exp(theta1(:,pnr,4)),'g-'), hold

pnr = 2;
% figure(2)
% % %xi = (-xiref + xiref*exp(theta1(:,pnr,1))/(1 + exp(theta1(:,pnr,1))));
% plot((-xiref + xiref*exp(theta1(:,pnr,1))./(1 + exp(theta1(:,pnr,1)))),'b-'), hold,
% plot((-xiref + xiref*exp(theta1(:,pnr,2))./(1 + exp(theta1(:,pnr,2)))),'k-')
% plot((-xiref + xiref*exp(theta1(:,pnr,3))./(1 + exp(theta1(:,pnr,3)))),'r-')
% plot((-0 + xiref*exp(theta1(:,pnr,4))./(1 + exp(theta1(:,pnr,4)))),'g-'), hold
psi2 = [theta1(:,pnr,1); theta1(:,pnr,2); theta1(:,pnr,3); theta1(:,pnr,4)];
figure(3), hist(psi2,100)

pnr = 2;
xi2a = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,1))./(1 + exp(theta1((burnin + 1):Lb,pnr,1)));
xi2b = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,2))./(1 + exp(theta1((burnin + 1):Lb,pnr,2)));
xi2c = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,3))./(1 + exp(theta1((burnin + 1):Lb,pnr,3)));
xi2d = -xiref + xiref*exp(theta1((burnin + 1):Lb,pnr,4))./(1 + exp(theta1((burnin + 1):Lb,pnr,4)));
xi2mcmc = [xi2a; xi2b; xi2c; xi2d];
figure(4), hist(xi2mcmc,100)
xihat = median(xi2mcmc);
xiint = prctile(xi2mcmc,[2.5 97.5]);

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
% The spatial dimensions of the rectangular area.
min_lon = min(longb) - 0.10;
max_lon = max(longb) + 0.05;
minlong = -23.0;
maxlong = -18.95;

mu = 4.5;
[sort_long,indx_long] = sort(longb);

figure(1), plot([minlong maxlong],ones(1,2)*(mu+(sigmahat/(-xihat))),'k-')
hold
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.05)^(-xihat) - 1)/(xihat))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.25)^(-xihat) - 1)/(xihat))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.5)^(-xihat) - 1)/(xihat))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.75)^(-xihat) - 1)/(xihat))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+(sigmahat*((1 - 0.95)^(-xihat) - 1)/(xihat))),'k--')
plot(sort(longb),mu + y(indx_long),'k.','MarkerSize',10)
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
axis([minlong maxlong (mu-0.05) (mu + 3.5)])
hold off
box on
print(1, '-dpdf', 'gp_const_magnitude_vs_longitude.pdf')

% Exponential distribution, constant sigma
minlong = -23.0;
maxlong = -18.95;
mu = 4.5;
[sort_long,indx_long] = sort(longb);

figure(2), hold
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.05)))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.25)))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.5)))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.75)))),'k--')
plot([minlong maxlong],ones(1,2)*(mu+((1/betahat)*(-log(1 - 0.95)))),'k--')
plot(sort(longb),mu + y(indx_long),'k.','MarkerSize',10)
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
axis([minlong maxlong (mu - 0.05) (mu + 3.5)])
hold off
box on
print(2, '-dpdf', 'expon_const_magnitude_vs_longitude.pdf')

grmat1a = reshape(theta1((burnin + 1):Lb,1,:),L,m);
[upper1a,neff1a,lag1_corr1a,acpt_rate1a] = gpar(grmat1a)

grmat1b = reshape(theta1((burnin + 1):Lb,2,:),L,m);
[upper1b,neff1b,lag1_corr1b,acpt_rate1b] = gpar(grmat1b)

% The updated Gutenberg-Richer plot - based on GP
mu = 4.49;
zt = mu + sigmahat*(0:0.001:(1/(-xihat)-0.001));
fzt = n*exp((1/(-xihat))*log(1+xihat*(zt - mu)/sigmahat));
indx = 1:n;
Nrv = n*(1 - indx/(n + 1));
alp2 = n:(-1):1;
bet2 = n - alp2 + 1;
mean2 = 1 - indx/(n + 1);
lower2 = betainv(0.025,alp2,bet2);
upper2 = betainv(0.975,alp2,bet2);
qtile1 = mu + sigmahat*(lower2.^(-xihat) - 1)./xihat;
qtile2 = mu + sigmahat*(mean2.^(-xihat) - 1)./xihat;
qtile3 = mu + sigmahat*(upper2.^(-xihat) - 1)./xihat;
ysc = y;
ysco = sort(ysc);
figure(16), semilogy((mu + ysco),Nrv,'k.')
hold
semilogy(zt,fzt,'k-')
semilogy([(mu + sigmahat/(-xihat)) (mu + sigmahat/(-xihat))],[0.001 1000],'k-.')
semilogy(qtile1,Nrv,'k--')
%semilogy(qtile2,Nrv,'g-')
semilogy(qtile3,Nrv,'k--')
hold
set(gcf,'PaperUnits','centimeters')
set(gcf,'papersize',[15,10])
set(gcf,'paperposition',[0,0,15,10])
set(gcf,'units','centimeters')
set(gcf,'Position',[0,0,15,10])
xlabel('Magnitude','interpreter','Latex','fontsize',17)
ylabel('$\log(N)$','interpreter','Latex','fontsize',17)
set(gca,'fontsize',13)
axs=axis
axis([mu (mu + 3.52) (0.74) (n+40)])
hold off
box on
print(16, '-dpdf', 'gp_constant_gut_rich_updated.pdf')

% The updated Gutenberg-Richer plot
mu = 4.49;
zt = mu + sigmahat*(0:0.001:(1/(0.1)-0.001));
fzt = n*exp(-betahat*(zt - mu));
indx = 1:n;
Nrv = n*(1 - indx/(n + 1));
alp2 = n:(-1):1;
bet2 = n - alp2 + 1;
mean2 = 1 - indx/(n + 1);
lower2 = betainv(0.025,alp2,bet2);
upper2 = betainv(0.975,alp2,bet2);
qtile1 = mu - (1/betahat)*log(lower2);
qtile2 = mu - (1/betahat)*log(mean2);
qtile3 = mu - (1/betahat)*log(upper2);
ysc = y;
ysco = sort(ysc);
figure(17), semilogy((mu + ysco),Nrv,'k.')
hold
semilogy(zt,fzt,'k-')
%semilogy([(mu + sigmahat/(-xihat)) (mu + sigmahat/(-xihat))],[0.001 1000],'k-.')
semilogy(qtile1,Nrv,'k--')
%semilogy(qtile2,Nrv,'g-')
semilogy(qtile3,Nrv,'k--')
hold
set(gcf,'PaperUnits','centimeters')
set(gcf,'papersize',[15,10])
set(gcf,'paperposition',[0,0,15,10])
set(gcf,'units','centimeters')
set(gcf,'Position',[0,0,15,10])
xlabel('Magnitude','interpreter','Latex','fontsize',17)
ylabel('$\log(N)$','interpreter','Latex','fontsize',17)
set(gca,'fontsize',13)
axs=axis
axis([mu (mu + 3.52) (0.74) (n+40)])
hold off
box on
print(17, '-dpdf', 'expon_constant_gut_rich_updated.pdf')

save('constant_gp_exp') 
