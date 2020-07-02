clc
clear all


%% From q-Method
% BbarN =
% 
%     0.4159   -0.8549    0.3101
%    -0.8338   -0.4946   -0.2453
%     0.3631   -0.1565   -0.9185
% q =
% 
%    -0.0264
%     0.8410
%    -0.5020
%     0.2001


%% Q.5 Find BbarN with QUEST, knowing v1 is the more accurate senso

% Attitude sensors in body frame
v1B = [0.8273, 0.5541, -0.0920]';
v2B = [-0.8285, 0.5522, -0.0955]';

% Attitude sensors in intertial frame
v1N = [-0.1517, -0.9669, 0.2050]';
v2N = [-0.8393, 0.4494, -0.3044]';

% Finding K matrix
w1 = 1;
w2 = 1;

I = eye(3);
B = w1 * v1B * v1N' + w2 * v2B * v2N';
S = B + B';
Z = [B(2,3)-B(3,2), B(3,1)-B(1,3), B(1,2)-B(2,1)]';
sigma = trace(B);

K = [sigma Z';
    Z S-sigma*I];

% Finding the characteristic equation
syms s;
f = matlabFunction(vpa(det(K - s * eye(4)), 5));
f_diff = matlabFunction(vpa(diff(det(K - s * eye(4))),5));

%Iteration by Newton-Raphson
n = 5; % Number of desired iterations
lambda = w1 + w2;
res =  zeros(n,1);

for i = 1:n
    res(i+1) = f(lambda);
    lambda = lambda - (f(lambda)/f_diff(lambda));
end


% Finding Classic Rodrigues Parameter (CRP)
disp('CRP (Not normalized)')
q = inv((lambda + sigma) * I - S) * Z

% Finding the Quartenion set
disp('Quaternion')
betha = 1/sqrt(1 + q'*q) * [1; q]

% Find BbarN
e = [betha(2), betha(3), betha(4)]';
e_tilde = [0, -e(3), e(2);
          e(3), 0, -e(1);
          -e(2), e(1), 0];
BbarN = (betha(1)^2 - e'*e)*I + 2*e*e' - 2*betha(1)*e_tilde