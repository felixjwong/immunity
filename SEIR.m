%% Simulates well-mixed SEIR model predictions.

%clear all;

AT = [];
BT = [];
SS = [];

for FF = 1:11
f = (FF-1)*0.001;

%% Parameters
beta = 0.02; % units 1/time, transmission rate
gamma = 0.0667; % units 1/time, recovery rate, transmission rate/recovery rate is basic reproduction number R0
delta = 1/5; % units 1/time, incubation rate
d = 10;

% initial S/E/I/R fractions
n = 10000;
E(1) = 0;
I(1) = n/100;
S(1) = n-I(1);
R(1) = 0;
totalI(1) = I(1); % total infected over time

% number of timesteps
dt = 1;
tend = 700;

reference_beta = beta;
pruned_beta_ss = beta;
for t=2:tend
    beta = pruned_beta_ss; 
    S(t) = S(t-1) - (beta*I(t-1)*S(t-1))*dt*d/n;
    E(t) = E(t-1) + beta*I(t-1)*S(t-1)*dt*d/n - delta*E(t-1)*dt;
    I(t) = I(t-1) + delta*E(t-1)*dt - gamma*I(t-1)*dt;
    if t == 70 %|| t == 220 || t == 370 || t == 520 || t == 670 || t == 820
        I(t) = I(t) + S(t-1)*f;
        S(t) = S(t) - S(t-1)*f;
        totalI(t) = totalI(t-1) + delta*E(t-1)*dt + S(t-1)*f;
    else
        totalI(t) = totalI(t-1) + delta*E(t-1)*dt;
    end
    
    R(t) = 1 - I(t) - S(t) - E(t); 
end

AT(FF)=(totalI(700)-reftotal(700))-(totalI(69)-reftotal(69));
BT(FF)=(totalI(70)-reftotal(70))-(totalI(69)-reftotal(69));
end

%figure;
%hold on; plot([0:0.1:1],AT/n,'k');plot([0:0.1:1],BT/n,'b')
%ylim([0 1])

%figure;plot(diff(totalI/n,1))
%xlim([0 700])

% figure; 
% hold on
 %plot(1+[1:length(totalI)]*dt,totalI/n,'--r','LineWidth',5)
% box on;
% xlabel('Time'); ylabel('Fraction')
% legend('Total Infected')
% ylim([0 1.2])
% xlim([0 tend*dt])