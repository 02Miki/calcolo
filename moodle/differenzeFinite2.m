%% 1 !!! controllare
clc
clear
close all


f = @(x) x+cos(x);

% SECONDA CONDIZIONE NON è 2/3, PERCHé è MU*DU/DX = 2, non du/dx
BC = [-1, 2];

a = 0;
b = pi;

mu = 3;

intervalli = 5;

h = (b-a)/intervalli;

x = (a+h:h:b)';
N = length(x);

uni = ones(N,1);

d = mu/h^2*2*uni;
d(end) = d(end)/2;

d1 = -mu/h^2*uni;

A = spdiags([d, d1, d1], [0 1 -1], N, N);

b = f(x);

b(1) = b(1) + mu/h^2 * BC(1);
b(end) = (b(end) + 2/h * BC(2))/2;

u = [BC(1); A\b];
max(u)
% 3.919 invece di 3.8729


%% 2

clc
clear
close all

f = @(x) exp(sin(x));

BC = [1, -2];

mu = @(x) exp(sin(x));
a = 0;
b = pi;
intervalli = 1000;
h = (b-a)/intervalli;

xMu = (a:h:b)';
% tolgo b perché so già quanto vale in quel punto (dirichlet)
x = a:h:b-h;

muVec = mu((xMu(1:end-1)+xMu(2:end))/2);

d =  1/h^2*[muVec(1); (muVec(1:end-1) + muVec(2:end))];

d1 = 1/h^2*(-muVec(1:end-1));

% spdiags vuole vettori colonna
A = spdiags([d, [0; d1], [d1; 0]], [0 1 -1], intervalli, intervalli);

terminiNoti = f(x)';

terminiNoti(1) = terminiNoti(1)/2 - BC(1)/h;
terminiNoti(end) = terminiNoti(end) + (muVec(end))*BC(2)/h^2;

u = [A\terminiNoti; BC(2)];

max(u)
% plot([x,b], u)


%% 3

clc
clear
close all


mu_f = @(x) x+1;

f = @(x) 2.*(x+1).*sin(x);

a = 0;
b = pi;
BC = [0 0];
N = 1000;

h = (b-a)/N;
x = (a:h:b)';

xM = (x(1:end-1) + x(2:end))/2; 
mu = mu_f(xM);


d = 1/h^2.*(mu(1:end-1) + mu(2:end));

d1 = -1/h^2.*mu(2:end-1);

A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N-1, N-1);

b = f(x(2:end-1));

b(1) = b(1) + mu(1)/h^2 * BC(1);
b(end) = b(end) - mu(end)/h^2 * BC(2);

u = [BC(1); A\b; BC(2)];
max(u)

% 2.1054, ok









