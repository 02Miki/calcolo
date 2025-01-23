%% 1
clc
clear
close all


A = [3.75 -1.25 -2.25 0.75; -1.25 3.75 0.75 -2.25; -2.25 0.75 3.75 -1.25; 0.75 -2.25 -1.25 3.75];

[v,d,f ] = eigs(A)


%% 2
clc
clear
close all

f = @(x) 1./(x.^2+1);
x = linspace(-1, 1,20+1)
trapz(x, f(x))


%% 3
clc
clear
close all

N = 2500;

A = spdiags(ones(N,1)*[2, -1 -1 -1 -1], [0 -500 -1 1 500], N,N);

% gauss sidel
% varia sia d che f (upper)

D = diag(diag(A));
F = triu(A) - D;
E = A-D-F;
b = ones(N,1);

x = zeros(N,1);
for k = 1:68
    x = (b-E.*x)\(F+D);
        k
    if norm(A*b-x)/norm(b) <= 10^-4
        k
        break
    end
end


%% 4
clc
clear
close all

f = @(x) 1./(x.^4+1);

a = -4;
b = 6;
N = 6+1;
x = linspace(a, b, N)

fit = polyfit(x, f(x), 6);

y = polyval(fit, linspace(a, b))

figure
hold on
plot(linspace(a, b), y)

plot(linspace(a,b), f(linspace(a,b)))

%% 5

clc
clear
close all

f = @(x) 1/3*x^3 + 3*x + 2*sin(x) + pi;

a = -1;
b = 0;

for k = 1:12
    xM = (a+b)/2
    if f(a)*f(xM) < 0
        b = xM;
    elseif f(b)*f(xM) <0
        a = xM;
    end
end



%% 6
clc
clear
close all

phi = @(x) log(x^2+1) +1;

x0 = 0;
for k = 1:1000
    x = phi(x0);

    if abs(x-x0) <= 10^-6
        k
        break
    end
    x0 = x;


end


%% 7 


clc
clear
close all

x = [0 0.1 0.12 0.24 0.31 0.35 0.41 0.44 0.56 0.67 0.73 0.85 0.92 1]';

N = length(x) - 2;

BC = [0 0];

mu_F = @(x) cos(x);

f = @(x) x+1;

h = diff(x);

xMedi = (x(2:end) + x(1:end-1))/2;

muVal = mu_F(xMedi);

d = muVal(1:end-1)./h(1:end-1) + muVal(2:end)./h(2:end);

d1 = -muVal(2:end-1)./h(2:end-1);

A = spdiags([d [0; d1] [d1;0]], [0 1 -1], N, N);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

% inutile perché omogeneo
b(1) = b(1) + muVal(1)/h(1) * BC(1);
b(end) = b(end) + muVal(end)/h(end) * BC(2);


u = [BC(1); A\b; BC(2)];


max(u)

% 0.2381, ok








