clc
clear all

%% From TRIAD
% BbarN =
% 
%     0.4156   -0.8551    0.3100
%    -0.8339   -0.4943   -0.2455
%     0.3631   -0.1566   -0.9185


%% Q.6 Find BbarN with q-method, knowing v1 is the more accurate sensor

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

% Find the optimal attitude (Max eigenvalue/eigenvector)
[V, v] = eig(K);
v = diag(v);
index_v = find(v==max(v));

q = V(:,index_v)

% Find BbarN
e = [q(2), q(3), q(4)]';
e_tilde = [0, -e(3), e(2);
          e(3), 0, -e(1);
          -e(2), e(1), 0];
BbarN = (q(1)^2 - e'*e)*I + 2*e*e' - 2*q(1)*e_tilde


