%% 1

x = 0:8;
y = [4.9 4.3 7.1 3.4 2.9 2.1 3.5 7.3 2.3];

splineVincolata = spline(x,[-0.5 y -2], 2.5);

splineVincolata
% 5.5823

%% 2
clc
clear
close all


N = 250000;
uni = ones(N, 1);

A = spdiags([4 * uni, -1*uni, -1*uni, -1*uni, -1*uni], [0, -500, -1, 1, 500], N, N);
b = uni;
x = A\b;
x(561)

% 497.2

%% 3
clc
clear 
close all


n = 5;
b = 1;
a = 0;
f = @(x) exp(-1./(x+1));
k = 0:n;
x = cos((2*k+1)/(2*(n+1))*pi);

% bisogna riscalare e traslare i nodi di chebychev

z = (b-a)/2.*x+(a+b)/2;

% in potenza discendente
fit = polyfit(z, f(z), n)

% -0.1861

%% 4
clc
clear
close all


n = 1:20;
A = zeros(20, 20);
b = ones(20, 1);
for r =n
   for c = n
        A(r, c) = 1+exp(-abs(r-c));
   end
   if rem(r, 2) == 0
       b(r,1) = 0;
   else
       b(r,1) = 1;
   end
end

x = A\b;
norm(x, Inf)

% 1.1034



%% 5

clc
clear
close all

f = @(x) 1./(x.^4+1);

n = 6;

a = -4;
b = 6;
x = linspace(a, b, n+1);
fit = polyfit(x, f(x), n);
val = polyval(fit, linspace(a, b, 500));

figure
hold on


plot(x, f(x), "b");
plot(linspace(a, b, 500), val)

% 5

%% 6
% RIFAI

clc
clear
close all


N = 2500;

uni = ones(N, 1);

A = spdiags([2*uni, -uni -uni -uni -uni], [0, -500, -1, 1, 500], N, N);


b = uni;

x = zeros(N, 1);

k = 0;
M_GS = tril(A);
N_GS = A - M_GS;
rho_GS = myRho(M_GS,N_GS)
while true

for n = 1:size(A, 1)
    % a = A(n,n)
    % laT = (b(n,:) - A(n,n+1:end) * x(n+1:end) - A(n, 1:n-1) * x(1:n-1))
    % laB= b(n,:)
    % dx = A(n,n+1:end) * x(n+1:end) - A(n, 1:n-1) * x(1:n-1)
    % soloX = x(n)
    x(n) = A(n,n)\(b(n,:) - A(n,n+1:end) * x(n+1:end) - A(n, 1:n-1) * x(1:n-1));
    % disp("=====")

end


    if norm(A*x-b)/norm(b) < 10^-4
        break
    end
    k = k+1;
    if k == 100
        break
    end
end



% FATTO A MODO MIO, NON CONVERGE

%% 6 FATTO BENE



%% 7

A = [3.75 -1.25 -2.25 0.75; -1.25 3.75 0.75 -2.25; -2.25 0.75 3.75 -1.25; 0.75 -2.25 -1.25 3.75];

[e, V] = eig(A)

% 3

%% 8

clc
clear
close all



spazio = linspace(0, 2*pi, 21)';
f = @(x)  x+sin(x).*cos(x)+0.4.*(1-sin(x).^2);

A = [ones(length(spazio), 1), spazio, cos(spazio), sin(spazio)];

AT = A';
y = AT*f(spazio);
V = AT*A;

c = V\y


% !!! viene uguale facendo A\f(spazio)














