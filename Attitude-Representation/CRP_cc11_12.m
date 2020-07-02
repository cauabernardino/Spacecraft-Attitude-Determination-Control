clc
clear all


%% Q.1 Equivalent DCM value
q = [0.1, 0.2, 0.3]';

qtq = q'*q;
I = eye(3);
q_tilde = [0, -q(3), q(2);
           q(3), 0, -q(1);
           -q(2), q(1), 0];

C = (1/(1+qtq))*((1-qtq)*I + 2*q*q' - 2*q_tilde)


%% Q.2 Equivalent CRP

D = [0.333333, -0.666667, 0.666667;
    0.871795, 0.487179, 0.0512821;
    -0.358974, 0.564103, 0.74359];

zeta = sqrt(trace(D)+1);

q2 = (1/zeta^2)*[D(2,3)-D(3,2);
                D(3,1)-D(1,3);
                D(1,2)-D(2,1)]
            
%q2_tilde = (D'-D)/zeta^2
 
 
%% Q.4 B/F in terms of CRPs

FN = [0.1, 0.2, 0.3]';
BN = [-0.3, 0.3, 0.1]';
NF = -1*FN;

% BN to DCM
BNTBN = BN'*BN;
BN_tilde = [0, -BN(3), BN(2);
           BN(3), 0, -BN(1);
           -BN(2), BN(1), 0];

BNmatrix = (1/(1+BNTBN))*((1-BNTBN)*I + 2*BN*BN' - 2*BN_tilde);

% NF to DCM
NFTNF = NF'*NF;
NF_tilde = [0, -NF(3), NF(2);
           NF(3), 0, -NF(1);
           -NF(2), NF(1), 0];
       
NFmatrix = (1/(1+NFTNF))*((1-NFTNF)*I + 2*NF*NF' - 2*NF_tilde);

BF = BNmatrix*NFmatrix;

% CRP from BF

zeta2 = sqrt(trace(BF)+1);

q4 = (1/zeta2^2)*[BF(2,3)-BF(3,2);
                BF(3,1)-BF(1,3);
                BF(1,2)-BF(2,1)]
       
  
