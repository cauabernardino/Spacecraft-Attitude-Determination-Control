clc
clear all

C = [0.925417, 0.336824, 0.173648;
     0.0296956, -0.521281, 0.852869;
     0.377786, -0.784102,-0.492404];

phi = acos(0.5 * (C(1,1) + C(2,2) + C(3,3) - 1))

e = 1/(2 * sin(phi)) * [C(2,3) - C(3,2);
                        C(3,1) - C(1,3);
                        C(1,2) - C(2,1)]
                    
PHI = rad2deg(phi)
