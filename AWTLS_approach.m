%% Initializing the parameters for our simulation

Q_initial = 100;            % true initial capacity [Ah]
I_max = 5*Q_initial;        % must be able to measure current up to +/- maxI
I_sensor_bit = 2^10;        % 10-bit precision on current sensor (number of combinations)
slope = -0.01;              % rate of change of capacity
Q_nom = .99*Q_initial;       % nominal capacity
SoC_max = 0.8;                 % max change in SOC between capacity estimates
SoC_min = -SoC_max;               % min change in SOC (signed)
mode = 0.5; sigma = 0.6;    % The mode and standard deviation of the measurement updates
socnoise = 0.01;            % standard deviation of SOC estimates
gamma = 0.98;               % forgetting factor (Range Between 0.9 -0.99 depending on the use case)

%% Computing the simulated data sets
n = 1000; % number of data points to use in computations/plots
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
figure(1);
clf
plot(abs(y)/abs(x));
title("Using regular regression");
xlabel('Cycles');
ylabel('Capacity[Ah]');

%% kick off the algorithm with the measured parameters and the variance of the measurements (SD^2)
SigmaY = del_y.^2;
SigmaX = del_x.^2;
measX = x;
measY = y;

%% AWLTS Start
Qhat = zeros(length(measX),4); SigmaQ = Qhat; Fit = Qhat; 
K = sqrt(SigmaX(1)/SigmaY(1));  
% Initialize some variables used for the recursive methods.
C1 = 0; C2 =0; C3 = 0; C4 = 0; C5 = 0; C6 = 0;
  if Q_nom ~= 0
    C1 = 1/(K^2*SigmaY(1)); C2 = K*Q_nom/(K^2*SigmaY(1)); 
    C3 = K^2*Q_nom^2/(K^2*SigmaY(1));
    C4 = 1/SigmaX(1); C5 = K*Q_nom/SigmaX(1); C6 = K^2*Q_nom^2/SigmaX(1);
  end
for iter = 1:length(measX)
    % Compute some variables used for the recursive methods
    C1 = gamma*C1 + measX(iter)^2/(K^2*SigmaY(iter));
    C2 = gamma*C2 + K*measX(iter)*measY(iter)/(K^2*SigmaY(iter));
    C3 = gamma*C3 + K^2*measY(iter)^2/(K^2*SigmaY(iter));
    C4 = gamma*C4 + measX(iter)^2/SigmaX(iter);
    C5 = gamma*C5 + K*measX(iter)*measY(iter)/SigmaX(iter);
    C6 = gamma*C6 + K^2*measY(iter)^2/SigmaX(iter);
    %AWLTS
    r = roots([C5 (-C1+2*C4-C6) (3*C2-3*C5) (C1-2*C3+C6) -C2]);
    r = r(r==conj(r)); % discard complex-conjugate roots
    r = r(r>0); % discard negative roots
    Jr = ((1./(r.^2+1).^2).*(r.^4*C4-2*C5*r.^3+(C1+C6)*r.^2-2*C2*r+C3))';
    J = min(Jr);
    Q = r(Jr==J); % keep Q that minimizes cost function
    H = (2/(Q^2+1)^4)*(-2*C5*Q^5+(3*C1-6*C4+3*C6)*Q^4+(-12*C2+16*C5)*Q^3 ...
          +(-8*C1+10*C3+6*C4-8*C6)*Q^2+(12*C2-6*C5)*Q+(C1-2*C3+C6));
    Qhat(iter) = Q/K;
    SigmaQ(iter) = 2/H/K^2; %the variance calculated for error bounds 
    Fit(iter) = gammainc(J/2,(2*iter-1)'/2,'upper');
end
Fit = real(Fit);
figure(2);
clf;
plot(Qhat+3*sqrt(SigmaQ),'linewidth',0.5);
hold on
plot(Qhat-3*sqrt(SigmaQ),'linewidth',0.5); 
hold on
plot(Qhat);
hold on
plot(1:n,Qreal,'k:','linewidth',1); 
title('AWTLS Method');
xlabel('Cycles');
ylabel('Capacity[Ah]');
figure(3); clf; hold on; plot(Fit);