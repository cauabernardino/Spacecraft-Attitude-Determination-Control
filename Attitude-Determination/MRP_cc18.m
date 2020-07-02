clc
clear all

%% Q.1 Map sigma to the corresponding DCM

sigma1 = [0.1, 0.2, 0.3]';

sigma_tilde = [0, -sigma1(3), sigma1(2);
               sigma1(3), 0, -sigma1(1);
               -sigma1(2), sigma1(1), 0];

I = eye(3);

disp('Q.1 The equivalent DCM is:')
C = I + (8*sigma_tilde^2 - 4*(1-sigma1'*sigma1)*sigma_tilde)/(1+sigma1'*sigma1)^2


%% Q.2 Map the DCM to the equivalent MRP set

D = [0.763314, 0.0946746, -0.639053;
    -0.568047, -0.372781, -0.733728;
    -0.307692, 0.923077, -0.230769];

zeta = sqrt(trace(D)+1);
disp('Q.2 The equivalent MRP set is:')
sigma2 = 1/(zeta*(zeta+2)) * [D(2,3)-D(3,2);
                              D(3,1)-D(1,3);
                              D(1,2)-D(2,1)]

