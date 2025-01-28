%% 1

clc
clear
close all


f = @(n) (-1).^n .* (n.^2 + 3.*n);


n = 0:1000;

sum(f(n))



%% 2

clc
clear
close all


phi = @(x) (exp(-x^2)*sin(5*x) + 0.3)/(-10)

x0 = -1;
for k = 1:1000
    
    x = phi(x0);
    
    if abs(x-x0) < 10^-3
        k
        break
    end
    x0 = x;
end


%% 3

clc
clear
close all


f = @(x) exp(-x.^2);

N = [4 8 16 32];

a = -1;
b = 1;

I = integral(f, a, b);
errore = [];

for k = N
    x = linspace(a, b, k+1);
    val = trapz(x, f(x));

    errore(end+1) = abs(val-I);

end

errore(2:end)./errore(1:end-1)



% si sarebbe potuto anche solo calcolare N = 32 e N = 16



%% 4

clc
clear
close all


n = 256;
alfa = 2;

uni = ones(n,1);
d = alfa * uni;

d1 = -1*uni;

A = spdiags([d, d1, d1], [0 1 -1], n, n);

b = sum(A, 2);

[x, flag, residuo, iter] = pcg(A,b)

% !! flag >= 1 non è andato a convergenza



%% 5

clc
clear
close all


f = @(t,z) [z(2); z(3); z(4); -4*z(3)-4*z(1) + t^2+1];
tspan = [1 6];
BC = [0; 2; 4; 8];


[t45,y45] = ode45(f,tspan,BC);

[t,y] = ode15s(f,tspan,BC);

passi45 = size(t45)
passi15s = size(t)

passi45 = t45(2:end)-t45(1:end-1)

passi15s = t(2:end)-t(1:end-1)

%% 6

clc
clear
close all

mu_f = @(x) 2+x;
L = 2;

f = @(x) 1 + 1* (x>L/2);

x = L*[0 0.02 0.12 0.19 0.28 0.33 0.41 0.42 0.48 0.5 0.53 0.56 0.61 0.69 0.74 0.78 0.83 0.88 0.94 0.99 1]';

h = diff(x);

xMedi = (x(2:end) + x(1:end-1))/2;

muVal = mu_f(xMedi);

d = muVal(2:end)./h(2:end) + muVal(1:end-1)./h(1:end-1);

d1 = -muVal(2:end-1)./h(2:end-1);

N = length(x)-2;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N, N);


b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

% bisognerebbe modificare b con dirichlet, ma è omogeneo quindi sommerei 0

u = [0; A\b; 0];

u(find(x==L/2))


%% 7

clc
clear
close all


f = @(x) exp(exp(x));

integral(f, 0, 1)

%% 8

clc
clear
close all


n = 2500;

uni = ones(n,1);

d = 5*uni;

A = spdiags(d, 0, n, n);

A(1:end, end) = 1;

A(end, 1:end) = 1;

A(end, end) = n;

% metodo jacobi, varia solo diagonale

% Dx+1 = b-(E+F)x

D = diag(diag(A));
E = tril(A) - D;
F = A-D-E;

b = uni;
x = zeros(n,1);

for k = 1:36
    
    x = D\(b-(E+F)*x);

    if (norm(A*x-b)/norm(b)) <= 10^-8
        k
        break
    end

end


