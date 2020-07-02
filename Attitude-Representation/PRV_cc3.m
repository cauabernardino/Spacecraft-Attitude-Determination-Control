clc
clear all

t1 = 20;
t2 = -10;
t3 = 120;

%% Question 2
% 3-2-1 Rotation
C = [cosd(t2)*cosd(t1), cosd(t2)*sind(t1), -sind(t2);
      sind(t3)*sind(t2)*cosd(t1)-cosd(t3)*sind(t1), sind(t3)*sind(t2)*sind(t1)+cosd(t3)*cosd(t1), sind(t3)*cosd(t2);
      cosd(t3)*sind(t2)*cosd(t1)+sind(t3)*sind(t1), cosd(t3)*sind(t2)*sind(t1)-sind(t3)*cosd(t1), cosd(t3)*cosd(t2)]

% Mapping to Principal Rotation Parameters
phi = acos(0.5 * (C(1,1) + C(2,2) + C(3,3) - 1))

e = 1/(2 * sin(phi)) * [C(2,3) - C(3,2);
                        C(3,1) - C(1,3);
                        C(1,2) - C(2,1)]                    
% PHI = rad2deg(phi)


%% Question 3
FB = [1 0 0; 0 0 1; 0 -1 0];
BN = FB;

FN = FB*BN

% Mapping to Principal Rotation Parameters
phi3 = acos(0.5 * (FN(1,1) + FN(2,2) + FN(3,3) - 1))

e = 1/(2 * sin(phi3)) * [FN(2,3) - FN(3,2);
                        FN(3,1) - FN(1,3);
                        FN(1,2) - FN(2,1)]
PHI3 = rad2deg(phi3) 
