clc
clear all

%% Q.2 Find sigma R/N

sigmaBN = [0.1, 0.2, 0.3]';
sigmaRB = [-0.1, 0.3, 0.1]';

% Transforming into DCM for ease purpose
I = eye(3);
sigmaBNtilde = [0, -sigmaBN(3), sigmaBN(2);
                sigmaBN(3), 0, -sigmaBN(1);
               -sigmaBN(2), sigmaBN(1), 0];

sigmaRBtilde = [0, -sigmaRB(3), sigmaRB(2);
                sigmaRB(3), 0, -sigmaRB(1);
               -sigmaRB(2), sigmaRB(1), 0];

BN = I + (8*sigmaBNtilde^2 - 4*(1-sigmaBN'*sigmaBN)*sigmaBNtilde)/(1+sigmaBN'*sigmaBN)^2;
RB = I + (8*sigmaRBtilde^2 - 4*(1-sigmaRB'*sigmaRB)*sigmaRBtilde)/(1+sigmaRB'*sigmaRB)^2;

% Find the DCM RN
RN = RB*BN

% Find sigmaRN
zeta = sqrt(trace(RN)+1);
sigmaRN = 1/(zeta*(zeta+2)) * [RN(2,3)-RN(3,2);
                              RN(3,1)-RN(1,3);
                              RN(1,2)-RN(2,1)]
                          
%% Find sigma B/R
% sigma BN is the same of last question

sigmaRN2 = [0.5, 0.3, 0.1]';
sigmaRN2tilde = [0, -sigmaRN2(3), sigmaRN2(2);
                sigmaRN2(3), 0, -sigmaRN2(1);
               -sigmaRN2(2), sigmaRN2(1), 0];
RN2 = I + (8*sigmaRN2tilde^2 - 4*(1-sigmaRN2'*sigmaRN2)*sigmaRN2tilde)/(1+sigmaRN2'*sigmaRN2)^2;

% Find the DCM BR
BR = BN*RN2'

% Find sigmaRN
zeta2 = sqrt(trace(BR)+1);
sigmaRN = 1/(zeta2*(zeta2+2)) * [BR(2,3)-BR(3,2);
                              BR(3,1)-BR(1,3);
                              BR(1,2)-BR(2,1)]


