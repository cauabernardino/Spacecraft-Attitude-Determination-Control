clc
clear all

tf = 60;
s = 0.01;
t = 0:s:tf;
N = length(t);

q0 = [0.4, 0.2, -0.1]';
I = eye(3);

B = @(q) 0.5 * (I + [0, -q(3), q(2);
                q(3), 0, -q(1);
                -q(2), q(1), 0] + q * q');
w = @(t) deg2rad(3*[sin(0.1*t), 0.01, cos(0.1*t)])';

qp = zeros(N, 3);

qp(1,1) = q0(1);
qp(1,2) = q0(2);
qp(1,3) = q0(3);


for k = 1:N-1
    
    qk = q0 + s * B(q0) * w(t(k));
    q0 = qk;
    qp(k+1,:) = q0;
 
end

R = sqrt(qp(4201,1)^2 + qp(4201,2)^2 + qp(4201,3)^2)


figure

plot(t,qp(:,1),'-','LineWidth',1.25);
hold on
plot(t,qp(:,2),'--','LineWidth',1.25);
plot(t,qp(:,3),'-.m','LineWidth',1.25);
hold off
legend('q1','q2', 'q3');
xlabel('Time')
ylabel('Angles(rad)')
grid;
