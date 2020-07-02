clc
clear all

%% Q.1 Finde BbarN with TRIAD method, knowing v1 is the more accurate sensor

% Attitude sensors in body frame
v1B = [0.8273, 0.5541, -0.0920]';
v2B = [-0.8285, 0.5522, -0.0955]';

% Attiude sensors in intertial frame
v1N = [-0.1517, -0.9669, 0.2050]';
v2N = [-0.8393, 0.4494, -0.3044]';


% Find BbarT
t1B = v1B/norm(v1B);
t2B = cross(v1B, v2B)/norm(cross(v1B, v2B));
t3B = cross(t1B, t2B);
 
BbarT = [t1B t2B t3B];


% Find NT
t1N = v1N/norm(v1N);
t2N = cross(v1N, v2N)/norm(cross(v1N, v2N));
t3N = cross(t1N, t2N);

NT = [t1N t2N t3N];

% Find BbarN
disp('Q.1')
BbarN = BbarT * NT'


%% Q.2 Find the estimation error between BbarN and BN in terms of principal
% rotation angles in units of degrees


BbarN2 = [0.969846, 0.17101, 0.173648;
         -0.200706, 0.96461, 0.17101;
         -0.138258, -0.200706, 0.969846];

BN2 = [0.963592, 0.187303, 0.190809;
      -0.223042, 0.956645, 0.187303;
      -0.147454, -0.223042, 0.963592];

disp('Q.2')
BbarB = BbarN2 * BN2'

phi = rad2deg(acos(0.5 * (BbarB(1,1) + BbarB(2,2) + BbarB(3,3) - 1)))