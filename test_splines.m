% A script to test the B-splines.

% The number of equally spaced interior knots.
kx = 5;
ky = 5;

% Delta x and Delta y.
dx = 1/(kx+1);
dy = 1/(ky+1);

% The order of the splines.
M = 4;

% Determine the number of functions.
Nx = kx + M;
Ny = ky + M;

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
x = 0.0:0.01:1.0;
y = 0.0:0.01:1.0;
lx = length(x);
ly = length(y);

% Compute the x-splines and the y-splines.
[X] = spline_functions(x,tau_x,dx,kx,M);
[Y] = spline_functions(y,tau_y,dy,ky,M);

figure(1), clf, hold
for s = 1:Nx
   plot(x,X(:,s))
end
title('B -spline functions')
xlabel('x')
ylabel('B_{i,4}(x)')
axis([0 1 0 1])
box on
figure(1),hold
%print -depsc f:\hydrology_meteorology\presenation_1st_nbbc_07\bspline_1dim

% Create a matrix with the tensor product of X and Y.
Z = kron(Y',X');

% figure(2), clf, hold
% for j = 1:Nx*Ny
%     z_temp = reshape(Z(j,:),[lx ly]);
%     figure(2),surf(x,y,z_temp),colorbar
% end
% hold

% figure(2), clf
j = 14
z_temp = reshape(Z(j,:),[lx ly]);
figure(2),surf(x,y,z_temp),colorbar
xlabel('x')
ylabel('y')
zlabel('B_{7,7}(x,y)')
view(3)
box off
grid off
%print -depsc f:\hydrology_meteorology\presenation_1st_nbbc_07\bspline_2dim


