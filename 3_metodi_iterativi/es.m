%% 3
clc
clear
close all

A=[3,0,4;
    7,4,2;
    -1,-1,-2]; b=[7;13;-4];

diagA = diag(diag(A));

% B = -D\(E+F)
raggioJ = mioRho(-diagA, A-diagA)
% -(D+E)\F
raggioGS = mioRho(-tril(A), A-tril(A))








