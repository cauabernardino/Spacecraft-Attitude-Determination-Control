%% Axis Attitude Control with MRPs

clc
clear all

%% Integration parameters
tf = 120;
h = 0.01;
t = 0:h:tf;
N = length(t);

%% Initial states

I = [100 0 0; 0 75 0; 0 0 80]; % kg*m^2
sigmaBN = [0.1, 0.2, -0.1]';
omegaBN = deg2rad([30, 10, -20]'); 
omegaRN = [0 0 0]'; % Reference is zero

delOmega = omegaBN - omegaRN;

K = 5; %Nm
P = 10 * eye(3); %Nms
L = [0 0 0]';

sigma = zeros(N,3);
omega = zeros(N,3);

% Initialization

sigma(1,:) = sigmaBN';
omega(1,:) = delOmega';

% Function: sigmaDot = 0.25 * B * omega
B = @(sigma) 0.25 * ((1 - sigma' * sigma) * I + ...
    2*[0, -sigma(3), sigma(2);
    sigma(3), 0, -sigma(1);
    -sigma(2), sigma(1), 0] + 2 * sigma * sigma');

omegaDot = @(sigma, omega) inv(I) * (-K*sigma - P*omega);

%% Integration
for k = 1:N-1
   
    sigmak = sigma(k,:)' + h * B(sigma(k,:)') * omega(k,:)';
    omegak = omega(k,:)' + h * omegaDot(sigmak, omega(k,:)');
    
    % Check the necessity for shadow set
    b = sigmak'*sigmak;
    if b >= 1
        sigmak = -sigmak/b;
    end
   
    sigma(k+1,:) = sigmak';
    omega(k+1,:) = omegak';
end

% Norm at 30 s.
R = norm(sigma(3002,:))


figure

plot(t,sigma(:,1),'-','LineWidth',1.25);
hold on
plot(t,sigma(:,2),'--','LineWidth',1.25);
plot(t,sigma(:,3),'-.m','LineWidth',1.25);
hold off
legend('sigma1','sigma2', 'sigma3');
xlabel('Tempo (s)');
ylabel('MRPs');
grid;