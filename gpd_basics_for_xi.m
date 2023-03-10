% Facts about the GPD (gen. Pareto dist.)

% Assume median = log(2)*mu = qhalf where mu = 0.5

% Assume the upper bound for this mu is Mmax - Mmin = 7.3 - 4.0 = qmax

% Then under the GPD xi = log(1 - qhalf/qmax)/(-log(0.5))

% Calcuations:
mu = 0.5;
Mmax = 8.5;
Mmin = 4.0
qhalf = log(2)*mu;
qmax = Mmax - Mmin;
xi = log(1 - qhalf/qmax)/(-log(0.5))
sigma = qmax*(-xi)
ref = (-xi)*3
