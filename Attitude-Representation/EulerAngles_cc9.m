clc
clear all

tf = 60;
s = 0.01;
t = 0:s:tf;

N = length(t);

x0 = deg2rad([40 30 80]');

B = @(theta)(1/cos(theta(2)))*[0 sin(theta(3)) cos(theta(3)); 
                               0 cos(theta(2))*cos(theta(3)) -sin(theta(3))*cos(theta(2)); 
                               cos(theta(2)) sin(theta(2))*sin(theta(3)) cos(theta(3))*sin(theta(2))];

w = @(t)deg2rad(20*[sin(0.1*t) 0.01 cos(0.1*t)]');

rates = zeros(N,3);
rates(1,1) = x0(1);
rates(1,2) = x0(2);
rates(1,3) = x0(3);

res = 0;

for k = 1:N-1
    
    %x0 = wrapToPi(x0);
    xk = x0 + s*B(x0)*w(t(k));
    x0 = xk;
    rates(k+1,:) = x0';
 
end

rate42 = rates(4201,:)
nrate = norm(rate42)

psi = rates(:,1);
tht = rates(:,2);
phi = rates(:,3);

figure

plot(t,psi,'-','LineWidth',1.25);
hold on
plot(t,tht,'--','LineWidth',1.25);
plot(t,phi,'-.m','LineWidth',1.25);
hold off
legend('psi','theta','phi');
xlabel('Time')
ylabel('Angles(rad)')
grid;
