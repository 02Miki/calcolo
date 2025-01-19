%% 1

N = 250000;

uni = ones(N, 1);
A = spdiags(uni*[4, -1 -1 -1 -1], [0, -500, -1 1 500], N, N);

x = A\uni;
x(561)

%% 2

clc
clear
close all


f = @(x) exp(-x.^2);

a = -1;
b = 1;

N = 120;

h = (b-a)/N;
x = a:h:b;

simp = h/6 * (f(x(1)) + 2 * sum(f(x(2:end-1))) + 4 * sum(f(x(1:end-1) + h/2)) + f(x(end))     )

%% 3

clc
clear
close all


x = 0:8;
y = [4.9 4.3 7.1 3.4 2.9 2.1 3.5 7.3 2.3];

sp = spline(x, [-0.5 y -2], 2.5)



%% 4
clc
clear
close all

nodi = 21;


a = 0;
b = 2*pi;

% due modi per trovare le x
% h = (1+1)/(nodi-1)
% x = (a+b)/2 + (b-a)/2 * (-1:h:1)
x = linspace(a,b,nodi);

f = @(x) x+sin(x).*cos(x)+0.4.*(1-sin(x).^2);

y = f(x)';
A = [ones(nodi, 1), x', cos(x)', sin(x)'];

AT = A';

% equivalenti su matlab, su carta bisogna fare la prima
(AT*A)\(AT*y)

A\y


%% 5

clc
clear
close all

f = @(x) exp(-x.^2).*sin(5.*x) + 0.3 + x;

x = linspace(-5, 5)

figure
hold on
plot(x, f(x))
plot(zeros(100,1), x)
plot(x, zeros(100,1))
% la radice è negativa

x = -1;
df = @(x) exp(-x.^2).*sin(5.*x)*(-2*x) + exp(-x.^2).*cos(5.*x)*5+1;
for k = 1:20
    x = x - f(x)/df(x);
    if abs(f(x)) < 10^-6
        k
        break
    end
end


%% 6
clc
clear
close all

f = @(t,y) [sin(y(2)); -cos(t*y(1))]

[t,y] = ode45(f,[0,2],[1; 1])











