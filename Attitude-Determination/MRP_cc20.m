clc
clear all

tf = 60;
s = 0.01;
t = 0:s:tf;
N = length(t);

% Initial values
sigma0 = [0.4, 0.2, -0.1]';

I = eye(3);

% Function: sigmap = 0.25 * B * w
B = @(sigma) 0.25 * ((1 - sigma' * sigma) * I + ...
    2*[0, -sigma(3), sigma(2);
    sigma(3), 0, -sigma(1);
    -sigma(2), sigma(1), 0] + 2 * sigma * sigma');
w = @(t) deg2rad(20*[sin(0.1*t), 0.01, cos(0.1*t)])';

sigmap = zeros(N, 3);
sigmap(1,1) = sigma0(1);
sigmap(1,2) = sigma0(2);
sigmap(1,3) = sigma0(3);

for k = 1:N-1
    
    sigmak = sigma0 + s * B(sigma0) * w(t(k));
    sigma0 = sigmak;
    
    if norm(sigma0) >= 1
        sigma0 = -sigma0/(sigma0'*sigma0);
    end
    sigmap(k+1,:) = sigma0;
    
 
end

R = norm(sigmap(4201,:))


figure

plot(t,sigmap(:,1),'-','LineWidth',1.25);
hold on
plot(t,sigmap(:,2),'--','LineWidth',1.25);
plot(t,sigmap(:,3),'-.m','LineWidth',1.25);
hold off
legend('sigma1','sigma2', 'sigma3');
xlabel('Time')
ylabel('Angles(rad)')
grid;
