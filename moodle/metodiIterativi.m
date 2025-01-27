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

[x,f,resrel, iter] = pcg(A, b, [], [], L,L')

% iter = 4 se pcg(L), dovrebbe venire 3
% iter = 2 se pcg come sopra

%% Round 2

%% 1
clc
clear
close all


n = 256;
alfa = 2;

uni = ones(n,1);

d = alfa*uni;
d1 = -1*uni;


A = spdiags([d, d1, d1], [0 1 -1], n,n);

b = sum(A,2);


[x,flag,relres,iter] = pcg(A, b, 10^-2, 100)

%% 2

clc
clear
close all


% il metodo di gauss sidel è sicuramente convergente se rho(B_GS) < 1


alfa = 0.5;
A = [alfa, 1, 0; 1 4 -alfa; 0 -alfa, 2];

D = diag(diag(A));
E = tril(A) - D;
F = A-D-E;

B = -(E+D)\F;

max(abs(eigs(B)))




%% 3

clc
clear
close all

A = [5 1 0; 1 5 -3; 0 1 8];

b = [6 3 9]';

D = diag(diag(A));
E = tril(A)-D;
F = triu(A)-D;

B = -D\(E+F);

max(abs(eigs(B)))

%% 4 - !!
clc
clear
close all

n = 150;

uni = ones(n,1);

A = spdiags(uni*[2 -1 -1], [0 1 -1], n, n);




% B_gs = -(D+E)\F


D = diag(diag(A));
E = tril(A)-D;
F = triu(A)-D;

BGS = -(D+E)\F;
ris = max(abs(eigs(spdiags(BGS))))


clear
n = 150;

uni = ones(n,1);
B = spdiags(uni*[4 -1 -1], [0 1 -1], n, n);
D = diag(diag(B));
E = tril(B)-D;
F = triu(B)-D;

BGS = -(D+E)\F;
% !!! Serve spdiags altrimenti non va a convergenza
ris = max(abs(eigs(spdiags(BGS))))

%% 5

clc
clear
close all

n = 256;
alfa = 2;

uni = ones(n,1);

A = spdiags(uni*[alfa, 1, 1], [0 1 -1], n,n);

A(end,1) = 0.01;
A(1,end) = 0.01;

b = sum(A,2);

L = ichol(A);

[sx,sflag,srelres,siter] = pcg(A,b)
[x,flag,relres,iter] = pcg(A,b,[], [], L, L')











