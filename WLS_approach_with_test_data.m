clc
clear all

load Charge_Capacity.mat;
Q_eg = Q_test;
figure(1)
clf;
subplot(2,1,1);
plot(Q_eg);
%% Create a best fit plot for the dataset
n = length(Q_eg);
x = 1:length(Q_eg);
y = Q_eg;

% Get coefficients of a line fit through the data.
coefficients = polyfit(x, y, 1);
% Create a new x axis with exactly 1000 points (or whatever you want).
xFit = linspace(min(x), max(x), n);
% Get the estimated yFit value for each of those 1000 new x locations.
Q_eg = polyval(coefficients , xFit);
[dydx, yint] = polyfit(xFit,Q_eg,1);
hold on
plot(Q_eg);
title('Capacity degradation Data otained from HNEI')
xlabel('Cycles');
ylabel('Capacity[Ah]');
xlim([0 1103]);

%% Initializing the parameters for our simulation
Q_initial = dydx(2);            % true initial capacity [Ah]
I_max = 5*Q_initial;        % must be able to measure current up to +/- maxI
I_sensor_bit = 2^10;        % 10-bit precision on current sensor (number of combinations)
slope = dydx(1);              % rate of change of capacity
Q_nom = 1*Q_initial;       % nominal capacity
SoC_max = 0.8;                 % max change in SOC between capacity estimates
SoC_min = -SoC_max;               % min change in SOC (signed)
mode = 0.5; sigma = 0.6;    % The mode and standard deviation of the measurement updates
socnoise = 0.01;            % standard deviation of SOC estimates
gamma = 0.95;               % forgetting factor

%% Computing the simulated data sets
Qreal = (Q_initial+slope*(1:n))';   % simulated true capacity as it (possibly) changes over time
x = ((SoC_max-SoC_min)*rand(n,1)+SoC_min);   % simulated true change in SOC between measurements (z2-z1)
y = Qreal.*x;                       % true accumulated ampere hours over that interval

%% Add in some noise to both variables. (Generate Lognormal Random Number)
binsize = 2*I_max/I_sensor_bit; 
mu = log(mode)+sigma^2;
m = 3600*lognrnd(mu,sigma,n,1);
del_y = binsize*sqrt(m/12)/3600;   % the standard deviation of y (Current measurement noise)
del_x = socnoise*ones(n,1);        % the standard deviation of x (SoC estimation noise) 
%% Simulated Measured add the deviation (noise)
x = x + del_x.*randn(n,1);         % add noise to true change in SOC 
y = y + del_y.*randn(n,1);         % add noise to accumulated ampere hours
%% To show the effect of using simple linear regression
subplot(2,1,2);
plot(abs(y)/abs(x));
hold on
plot(1:n,Qreal,'k:','linewidth',1); 
title("Using regular regression");
xlabel('Cycles');
ylabel('Capacity[Ah]');
xlim([0 1103]);
%% kick off the algorithm with the measured parameters and the variance of the measurements (SD^2)
SigmaY = del_y.^2;
SigmaX = del_x.^2;
measX = x;
measY = y;
%% Wls start
Qhat = zeros(length(measX),4); SigmaQ = Qhat; Fit = Qhat; 
K = sqrt(SigmaX(1)/SigmaY(1));
c1 = 0; c2 =0; c3 = 0;
if Q_nom ~= 0
    c1 = 1/SigmaY(1); c2 = Q_nom/SigmaY(1); c3 = Q_nom^2/SigmaY(1);
end

for iter = 1:length(measX)
    % Compute some variables used for the recursive methods
    c1 = gamma*c1 + measX(iter)^2/SigmaY(iter);
    c2 = gamma*c2 + measX(iter)*measY(iter)/SigmaY(iter);
    c3 = gamma*c3 + measY(iter)^2/SigmaY(iter);
    
    % WLS
    Q = c2./c1; 
    Qhat(iter,1) = Q;
    H = 2*c1;   
    SigmaQ(iter,1) = 2/H;
    J = Q.^2.*c1 -2*Q.*c2 + c3;
    Fit(iter,1) = gammainc(J/2,(iter-1)/2,'upper');
end
Fit = real(Fit);
figure(2);
clf;
subplot(2,1,1);
plot(Qhat+3*sqrt(SigmaQ),'linewidth',0.5);
hold on
plot(Qhat-3*sqrt(SigmaQ),'linewidth',0.5); 
hold on
plot(Qhat,'LineWidth',1);
hold on
plot(1:n,Qreal,'k:','linewidth',1); 
title("Using WLS");
xlabel('Cycles');
ylabel('Capacity[Ah]');
xlim([0 1103]);
ylim([1 2.6])
subplot(2,1,2); 
hold on; 
plot(Fit,'LineWidth',2);
title('Goodness of fit');
xlim([0 1103]);