%% 1

clc
clear
close all


mu = 3;

f = @(x) x+cos(x);

% neumann a dx
BC = [-1 2];

intervalli = 1000;

a = 0;
b = pi;

h = (b-a)/intervalli;

% solo i nodi interni sono uguali
uni = ones(intervalli-1,1);

d = [mu/h^2 * 2 * uni; mu/h^2 *(1)];

d1 = mu/h^2 * (-1) * [uni; 1];

A = spdiags([d, d1, d1], [0 1 -1], intervalli, intervalli);

x = (a:h:b)';
b = f(x(2:end));

b(1) = b(1) + mu/h^2*BC(1);

b(end) = b(end)/2 + BC(2)/h;


u = [BC(1); A\b];

max(u)


%% 2

clc
clear
close all

a = 0;
b = pi;

f = @(x) exp(sin(x));

% neumann a sx
BC = [1 -2];

intervalli = 1000;

h = (b-a)/intervalli;

x = (a:h:b)';

mu_f = @(x) exp(sin(x));
xMedi = (x(1:end-1) + x(2:end))/2;
muVal = mu_f(xMedi);

d = muVal(1:end-1) + muVal(2:end);

% aggiungo extra neumann a sx
% d = [(mu_f(x(1)+h/2) + mu_f(x(1)-h/2))/2; d]./h^2;
d = [muVal(1); d]./h^2;
% non mi serve il valore finale visto che conosco già l'ultima u grazie a
% dirichlet
d1 = -muVal(1:end-1)./h^2;


A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], intervalli, intervalli);

b = f(x(1:end-1));

b(1) = b(1)/2 - BC(1)/h;

b(end) = b(end) + muVal(end) * BC(end)/h^2;


u = [A\b;BC(end)];
plot(x, u)
max(u)





%% 3
clc
clear
close all

a = 0;
b = pi;

f = @(x) 2.*(x+1).*sin(x);
mu_f = @(x) x+1;

BC = [0 0];

intervalli = 1000;

h = (b-a)/intervalli;
x = (a:h:b)';
xMedi = (x(1:end-1) + x(2:end))/2;

muVal = mu_f(xMedi);

d = 1/h^2*(muVal(1:end-1) + muVal(2:end));

d1 = -muVal(2:end-1)/h^2;

A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], intervalli-1, intervalli-1);

b = f(x(2:end-1));

max(A\b)


%% Seconda parte =======================================

%% 1

clc
clear
close all


mu = 1/2;

f = @(x) x-3*cos(x);

a = 0;
b = pi;

% sarebbe 0/mu, visto che nella traccia c'è u' = 0
% neumann a sinistra
BC = [0 1];

intervalli = 1000;
h = (b-a)/intervalli;
% -1 per prendere i nodi interni
uni = ones(intervalli-1, 1);


d = [mu/h^2; 2*mu/h^2 .* uni];

d1 = -mu/h^2*uni;

A = spdiags([d,[0;d1],[d1;0]], [0 1 -1], intervalli, intervalli);

x = (a:h:b)';
b = f(x(1:end-1));

b(1) = b(1) - BC(1)/h; 

b(end) = b(end) + BC(2)*mu/h^2;

max(A\b)


%% 2

clc
clear
close all

a = 0;
b = pi;


intervalli = 1000;

mu_f = @(x) 1+log(x+1);

f = @(x) sin(cos(x));

h = (b-a)/intervalli;
x = (a:h:b)';

xMedi = (x(1:end-1) + x(2:end))/2;

muVal = mu_f(xMedi);

d = (muVal(1:end-1) + muVal(2:end))./h^2;

d1 = -muVal(2:end-1)./h^2;

% -1 perché DD
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], intervalli-1, intervalli-1);

b = f(x(2:end-1));

max(A\b)


%% 3
clc
clear
close all

a = 0;
b = pi;

% neumann a dx
BC = [1 -1/2];

mu_f = @(x) exp(sin(x));

f = @(x) cos(cos(x) - sin(x));

intervalli = 5;
h = (b-a)/intervalli;
x = (a:h:b)';

xMedi = (x(1:end-1) + x(2:end))/2;

muVal = mu_f(xMedi);

d = muVal(1:end-1) + muVal(2:end);

d = [d; muVal(end)]./h^2;

d1 = -muVal(2:end)./h^2;

A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], intervalli, intervalli);

b = f(x(2:end));

b(1) = b(1) + BC(1)*muVal(1)/h^2;

b(end) = b(end)/2 + BC(2)/h;

max(A\b)

















