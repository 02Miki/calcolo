%% 1

clc
clear
close all

x = [0 1 2 3]';
y = [1 2 4 8]';


% a0*1+a1*e^x
A = [ones(length(x),1) exp(x)]

x = A\y

norm(A*x-y)
% !!!!!! chidere, viene diverso (dovrebbe venire 2.84e-1, viene 0.532830789423654)

%% 3
clc
clear
close all

x = [-4 6 16 23];
y = [5 3 4 9];

spline(x, [5 y 6], log(1.3))


