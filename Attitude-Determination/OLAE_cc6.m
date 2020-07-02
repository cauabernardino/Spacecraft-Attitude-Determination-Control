clc
clear all

%% Q.2 Find BbarN with OLAE method

% Attitude sensors in body frame
v1B = [0.8273, 0.5541, -0.0920]';
v2B = [-0.8285, 0.5522, -0.0955]';

% Attitude sensors in intertial frame
v1N = [-0.1517, -0.9669, 0.2050]';
v2N = [-0.8393, 0.4494, -0.3044]';


% Finding the parameters d, S and W
d1 = v1B - v1N;
d2 = v2B - v2N;
d = [d1; d2];

s1 = v1B + v1N;
s2 = v2B + v2N;

s1_tilde = [0, -s1(3), s1(2);
          s1(3), 0, -s1(1);
          -s1(2), s1(1), 0];

s2_tilde = [0, -s2(3), s2(2);
          s2(3), 0, -s2(1);
          -s2(2), s2(1), 0];
      
S = [s1_tilde; s2_tilde];

index = size(d);

W = eye(index(1));


% Finding qBar
disp('CRP (Not normalized)')
qBar = inv(S' * W * S) * S' * W * d


% Finding the Quartenion set
disp('Quaternion')
beta = 1/sqrt(1 + qBar'*qBar) * [1; qBar]

% Find BbarN
I = eye(3);
disp('BbarN (Rotation matrix)')
e = [beta(2), beta(3), beta(4)]';
e_tilde = [0, -e(3), e(2);
          e(3), 0, -e(1);
          -e(2), e(1), 0];
BbarN = (beta(1)^2 - e'*e)*I + 2*e*e' - 2*beta(1)*e_tilde
