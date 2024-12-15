%% 2

clc
clear
close all


A = [1 -9 1; -7 9 -4; -6 4 -7];

b = [3; 10; -7];

x = zeros(3,1);

D = diag(diag(A));
EF = A-D;

for k = 1:4
    x = D\(b-EF*x);
    

end

x(2)


%% 3
clc
clear
close all


N = 150;
uni = ones(N, 1);
A = spdiags([4*uni, -uni, -uni], [0, -1 1], N, N);


b = (1:150).^(1/3);

x = zeros(N,1);

D = diag(diag(A));
EF = A-D;
k = 0

while norm(A*x-b,1)/norm(b,1) >= 10^-7
    x = D\(b-EF*x);
    
    k = k+1;
end


k


%% 5
clc
clear
close all

N = 256;
a = 2;

uni = ones(N,1);
A= spdiags([a*uni, uni, uni], [0 1 -1], N, N);
A(1, N) = 0.01;
A(N,1) = 0.01;

b = sum(A, 2);

L = ichol(A);

[x,f,resrel, iter] = pcg(L, b)

% iter = 4, dovrebbe venire 3 !!!!!!CHIEDERE!!!!!!



