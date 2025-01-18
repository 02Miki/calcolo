%% 1
clc
clear
close all


f = @(x) (x-1).*sin(x)-2.*cos(x);

BC = [0, 0];

mu =1;
h = 0.1;
x = [0:h:1]';
% -2 perché non voglio gli estremi
N = length(x)-2;
uni = ones(N, 1);
d = mu/h^2*2*uni;
d1 = -mu/h^2*uni;
d_1 = d1;

A = spdiags([d, d1, d_1], [0 1 -1], N, N)

% xMedio = (x(1:end-1)+x(2:end))/2;

b = f(x(2:end-1));

b(1) = b(1) + 1/h^2*mu*BC(1);
b(end) = b(end) + 1/h^2*mu*BC(end);


figure
hold on
u = [BC(1); A\b; BC(end)]

plot(x, u)

uvera = @(x) 1/mu .* (x-1).*sin(x);

plot(x, uvera(x))
legend("calcolata", "vera")

%% 2 -- Uguale a 1 ma elementi finiti
clc
clear
close all

f = @(x) (x-1).*sin(x)-2.*cos(x);

BC = [0, 0];

mu = 1;

x = [0, 0.1, 0.15, 0.25, 0.45, 0.6, 0.7, 0.95, 1]';

h = x(2:end) - x(1:end-1);

d = mu * (1./h(1:end-1) + 1./h(2:end));

% uguali
d_1 = -mu./h(1+1:end-1);
d1 = -mu./h(1+1:end-1);

N = length(x)-2;

A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N, N);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));
% non necessario in questo caso

b(1) = b(1) + mu(1)/h(1)*BC(1);
b(end) = b(end) + mu(end)/h(end)*BC(2);


figure
hold on
u = [BC(1); A\b; BC(2)]

plot(x, u)

uvera = @(x) 1/mu .* (x-1).*sin(x);

plot(x, uvera(x))

%% 3
clc
clear
close all

f = @(x) 2.*exp(x).*sin(x);

mu = 1;
h = 0.01;
x = 0:h:1;
BC = [1, exp(1)*cos(1)];
% -2 perché non voglio gli estremi
N = length(x)-2;
% -1 perché non deve prendere la prima riga, che sarà diversa a causa di
% neumann
uni = ones(N,1);
% tutta la diagonale
d = mu/h^2*[1; 2*uni];


d1 = [-mu/h^2; mu/h^2*-uni];

d_1 = mu/h^2*-uni;

% 2 nella f perché:
% 1 lo calcolo prima
% end - 1 perché so già u quanto vale in end (BC(2))
b = [1/2 * (f(x(1))-2/h*BC(1)), f(x(2:end-1))];

b(end) = b(end) + mu/h^2*BC(2);

% NELLA DIAGONALE 1, SPDIAGS MANGIA SOPRA, IN QUELLA -1 SPDIAGS MANGIA SOTTO
% CHECK DOCUMENTAZIONE
% aggiungo elemento a d_1 per avere dimensioni uguali
A = spdiags([d, d1, [d_1; 0]], [0 1 -1], N+1, N+1);


u = [A\b'; BC(end)];


figure
hold on
plot(x, u)

uvera = @(x) 1/mu .* exp(x) .*cos(x);

plot(x, uvera(x))
legend("calcolata", "vera")

%% 3 elementi finiti

clc
clear
close all


f = @(x) 2.*exp(x).*sin(x);

mu = 1;
h = 0.1;
x = (0:h:1)';
% neumann all'inizio
BC = [1, exp(1)*cos(1)];

N = length(x)-1;
uni = ones(N, 1);
d = 2*mu/h*uni;
d(1) = d(1)/2;

d1 = -mu/h*uni;

Ao = spdiags([d, d1, d1], [0 1 -1], N, N)


bo = f(x(1:end-1))/2 * 2*h;
bo(end) = bo(end) + mu/h*BC(2);

bo(1) = (f(x(1)) * h/2)-BC(1);

u = [Ao\bo; BC(2)];


figure
hold on
plot(x, u)

uvera = @(x) 1/mu .* exp(x) .*cos(x);

plot(x, uvera(x))

