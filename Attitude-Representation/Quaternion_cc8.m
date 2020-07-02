clc
clear all

tf = 60;
s = 0.01;
t = 0:s:tf;
N = length(t);

q0 = [0.408248, 0., 0.408248, 0.816497]';
w = [sin(0.1*t), 0.01, cos(0.1*t)] * 20;


B = @(b) [b(1), -b(2), -b(3), -b(4);
          b(2), b(1), -b(4), b(3);
          b(3), b(4), b(1), -b(2);
          b(4), -b(3), b(2), b(1)];
w = @(t)deg2rad(20*[0, sin(0.1*t), 0.01, cos(0.1*t)])';

bp = zeros(N, 4);

bp(1,1) = q0(1);
bp(1,2) = q0(2);
bp(1,3) = q0(3);
bp(1,4) = q0(4);

for k = 1:N-1
    
    qk = q0 + s * 0.5 * B(q0) * w(t(k));
    q0 = qk;
    bp(k+1,:) = q0;
 
end

R = sqrt(bp(4201,2)^2 + bp(4201,3)^2 + bp(4201,4)^2)


figure

plot(t,bp(:,1),'-','LineWidth',1.25);
hold on
plot(t,bp(:,2),'--','LineWidth',1.25);
plot(t,bp(:,3),'-.m','LineWidth',1.25);
plot(t,bp(:,4),'-*k','LineWidth',1.25);
hold off
legend('q0','q1','q2', 'q3');
xlabel('Time')
ylabel('Angles(rad)')
grid;
