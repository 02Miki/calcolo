%% 1 esercizio con richiamo
clc
clear
close all

intervalli = 110;

BC = [0 -15*log(5)];

a = 0;
b = 1;

mu = 1;

h = (b-a)/intervalli;

uni = ones(intervalli+1, 1);

d = 2*mu/h.*uni;
d(1) = d(1)/2;
d(end) = d(end)/2;

d1 = -mu/h.*uni;

x = linspace(a, b, intervalli+1)';

A = spdiags([d, d1, d1], [0 1 -1], intervalli+1, intervalli+1);


gamma_f = @(x) 9.*x.^4.*log(5)^2 - 6.*x.*log(5);

xMedi = (x(2:end) + x(1:end-1))/2;

d = h/3 * (gamma_f(xMedi(1:end-1)) + gamma_f(xMedi(2:end)) );
d = [h/3*gamma_f(xMedi(1)); d; h/3*gamma_f(xMedi(end))];
d1 = h/6*gamma_f(xMedi(1:end));

A_richiamo = spdiags([d,  [0;d1], [d1;0]], [0 1 -1], intervalli+1, intervalli+1);

A_tot = A+A_richiamo;

b = zeros(length(x),1);

b(1) = b(1)-BC(1);
b(end) = b(end)+BC(end);

u = A_tot\b;

uex = @(x) 5.^(2-x.^3);

abs(u(end) - uex(1))


%% 1 test
clc
clear
close all
intervalli = 110;

BC = [0 -15*log(5)];

a = 0;
b = 1;

mu = 1;
gamma_f = @(x) 9.*x.^4.*log(5)^2 - 6.*x.*log(5);

N = intervalli+1
h = (b-a)/(N-1);
x = linspace(a,b,N)';
uni = ones(N,1);
d = 2*mu/h*uni;
d(1) = d(1)/2;
d(end) = d(end)/2;

d1 = -mu/h*uni;


A = spdiags([d, d1, d1], [0 1 -1], N, N);

xMedi = (x(2:end) + x(1:end-1))/2;

dgamma = 1/3 * (gamma_f(xMedi(1:end-1)) + gamma_f(xMedi(2:end)))*h;
dgamma = [gamma_f(xMedi(1))/3*h; dgamma; gamma_f(xMedi(end))*h/3];


dgamma1 = 1/6 * gamma_f(xMedi(1:end)) * h;

A_gamma = spdiags([dgamma, [0;dgamma1], [dgamma1;0]], [0 1 -1], N, N);

b = zeros(N, 1);

b(1) = b(1) * h/2 - BC(1);
b(end) = b(end) * h/2 + BC(2);



u = (A+A_gamma)\b;
figure
hold on

plot(x,u, "r")

u_esatta = @(x) 5.^(2-x.^3);

abs(u(end) - u_esatta(1))
plot(x, u_esatta(x))

%% 2
clc
clear
close all


f = @(x) 3.*cos(x)+ (1-x).*cos(x);

mu = 1;

a = 1;
b = 2;


x = [1 1.19 1.22 1.39 1.47 1.5 1.66 1.79 1.81 1.92 2]';

% neumann a dx
BC = [2 0];

h = diff(x);

d = mu./h(1:end-1) + mu./h(2:end);

d = [d; mu/h(end)];

% !!! fare attenzione a mettere il . !!!
d1 = -mu./h(2:end);

nodi = length(x)-1;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], nodi, nodi);

b = f(x(2:end-1))/2 .* ((h(1:end-1) + h(2:end)));

b = [b; f(x(end)) * h(end)/2 + BC(2)];

b(1) = b(1) + mu/h(1)*BC(1);

u = [BC(1); A\b];

u(find(x==1.47))

%% 3

clc
clear
close all

mu = 1;

BC = [0 0];


f = @(x) 3.*cos(x)+ (3-x.^3).*cos(x); 

x = [0 0.39 0.4 0.61 0.82 1.09 1.23 1.41 1.69 1.9 2]';

h = diff(x);

d = mu./h(1:end-1) + mu./h(2:end);

d1 = -mu./h(2:end-1);

nodi = length(x)-2;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], nodi, nodi);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

u = [BC(1); A\b; BC(2)];

u(find(x==0.82))

%% 4

clc
clear
close all


mu = 1;
a = 0;
b = 2;

BC = [0 0];

f = @(x) 2.*sin(x)+sqrt(x+2).*cos(x);

x = [0 0.27 0.53 0.65 0.88 1.07 1.26 1.59 1.63 1.92 2]';

h = diff(x);

d = mu./h(1:end-1) + mu./h(2:end);

d1 = -mu./h(2:end-1);

nodi = length(x)-2;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], nodi, nodi);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

u = [BC(1); A\b; BC(2)];

u(find(x==1.59))













































