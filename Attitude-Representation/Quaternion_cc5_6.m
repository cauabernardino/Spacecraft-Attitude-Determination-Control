clc 
clear all


%% Q.1
q = [0.235702, 0.471405, -0.471405, 0.707107];

q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

C = [q0^2+q1^2-q2^2-q3^2, 2*(q1*q2+q0*q3), 2*(q1*q3-q0*q2);
    2*(q1*q2-q0*q3), q0^2-q1^2+q2^2-q3^2, 2*(q2*q3+q0*q1);
    2*(q1*q3+q0*q2), 2*(q2*q3-q0*q1), q0^2-q1^2-q2^2+q3^2]


%% Q.3

B = [-0.529403, -0.467056, 0.708231;
    -0.474115, -0.529403, -0.703525;
    0.703525, -0.708231, 0.0588291];

b0 = 0.5 * sqrt(B(1,1)+B(2,2)+B(3,3)+1);
b1 = (B(2,3)-B(3,2))/(4*b0);
b2 = (B(3,1)-B(1,3))/(4*b0);
b3 = (B(1,2)-B(2,1))/(4*b0);

b = [b0 b1 b2 b3]

%% Q.4

t1 = 20;
t2 = 10;
t3 = -10;

%% Question 2
% 3-2-1 Rotation
D = [cosd(t2)*cosd(t1), cosd(t2)*sind(t1), -sind(t2);
      sind(t3)*sind(t2)*cosd(t1)-cosd(t3)*sind(t1), sind(t3)*sind(t2)*sind(t1)+cosd(t3)*cosd(t1), sind(t3)*cosd(t2);
      cosd(t3)*sind(t2)*cosd(t1)+sind(t3)*sind(t1), cosd(t3)*sind(t2)*sind(t1)-sind(t3)*cosd(t1), cosd(t3)*cosd(t2)]
  
d0 = 0.5 * sqrt(D(1,1)+D(2,2)+D(3,3)+1);
d1 = (D(2,3)-D(3,2))/(4*d0);
d2 = (D(3,1)-D(1,3))/(4*d0);
d3 = (D(1,2)-D(2,1))/(4*d0);

d = [d0 d1 d2 d3]
