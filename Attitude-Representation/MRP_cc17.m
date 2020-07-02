clc
clear all

%% Q.1 Find short or long rotation

sigma1 = [0.1, 0.2, 0.3];
norm(sigma1)


%% Q.2 Find short or long rotation

sigma2 = [1.2, -0.1, -0.001];
norm(sigma2)

%% Q.3 Find corresponding Shadow Set

sigma3 = [0.1, 0.2, 0.3];

sigma3_2 = dot(sigma3, sigma3)
sigma3_s = [-sigma3(1)/sigma3_2, -sigma3(2)/sigma3_2, -sigma3(3)/sigma3_2]
