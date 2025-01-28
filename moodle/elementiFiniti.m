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

%% round 2

%% 1

clc
clear
close all

mu = 1;
a = 0;
b = 1;

f = @(x) 3.*cos(x) + (3-x.^3).*sin(x);

x = [0 0.11 0.25 0.36 0.49 0.52 0.6 0.77 0.82 0.92 1]';

h = diff(x);

d = mu./h(1:end-1) + mu./h(2:end);

d1 = -mu./h(2:end-1);

N = length(x)-2;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N, N);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

% non serve modificare B visto che ho dirichlet omogeneo

u = [0; A\b; 0];

u(find(x==0.36))


%% 2

clc
clear
close all

f = @(x) 3.*cos(x) + (3-x.^3).*cos(x);

x = [1 1.25 1.43 1.62 1.98 2.07 2.35 2.54 2.67 2.82 3]';

h = diff(x);

mu = 1;

d = mu./h(1:end-1) + mu./h(2:end);

d1 = -mu./h(2:end-1);

N = length(x)-2;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N, N);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

u = [0; A\b; 0];

u(find(x==1.43))


%% 3
clc
clear
close all


mu = 1;

f = @(x) 2.*sin(x) + (1-x).*sin(x);

% neumann a dx
BC = [0 2];

x = [0 0.3 0.54 0.73 0.95 1.01 1.24 1.53 1.76 1.95 2]';

h = diff(x);

d = mu./h(1:end-1) + mu./h(2:end);

d = [d; mu./h(end)];

d1 = -mu./h(2:end);

N = length(x)-1;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N, N);


b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

b = [b; f(x(end))/2 * h(end) + BC(2)];

u = [BC(1); A\b];

u(find(x==1.01))


%% 4

clc
clear 
close all

BC = [1 -1];

gamma_f = @(x) ((1-2.*x).^2 - 2);

mu =1;

intervalli = 140;
a = 0;
b = 1;
x = linspace(a, b, intervalli+1)';

h = 1/intervalli;

d = 2*mu/h * ones(length(x)-2,1);

d = [mu./h; d; mu./h];

d1 = -mu./h*ones(length(x)-1, 1);

A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], intervalli+1, intervalli+1);

xMedi = (x(2:end) + x(1:end-1))/2;

d_gamma = h/3 * (gamma_f(xMedi(1:end-1)) + gamma_f(xMedi(2:end)));

d_gamma = [h/3 * gamma_f(xMedi(1)); d_gamma; h/3 * gamma_f(xMedi(end))];
d1_gamma = 1/6 *gamma_f(xMedi(1:end)) *h;

A_gamma = spdiags([d_gamma, [0;d1_gamma], [d1_gamma;0]], [0 1 -1], intervalli+1, intervalli+1);


A_tot = A+A_gamma;

b = zeros(intervalli+1, 1);

b(1) = -BC(1);
b(end) = +BC(end);
u = A_tot\b;

u_esatta = @(x) exp(x-x.^2);


figure 
hold on

plot(x, u, LineWidth=2)
plot(x, u_esatta(x), "o")



abs(u(end)-u_esatta(1))


%% tema esame

clc
clear
close all

L = 2;
x = L.*[0 0.02 0.12 0.19 0.28 0.33 0.41 0.42 0.48 0.53 0.56 0.61 0.69 0.74 0.78 0.83 0.88 0.94 0.99 1]';

f = @(x) 1+1.*(x>L/2);

mu_f = @(x) 2+x;
h = diff(x);

xMedi = (x(2:end) + x(1:end-1))/2;

muVal = mu_f(xMedi);

d = muVal(1:end-1)./h(1:end-1) + muVal(2:end)./h(2:end);

d1 = -muVal(2:end-1)./h(2:end-1);

N = length(x)-2;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N, N);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

u = [0; A\b; 0];

plot(x,u)
% più di 0.2466, meno di 0.24899 -> 0.247, OK


%% extra

%% 2

clc
clear
close all

f = @(x) 2.*sin(x) + (1-x).*sin(x);
x = [0 0.1 0.27 0.37 0.41 0.54 0.66 0.77 0.8 0.93 1]';

h = diff(x);

mu = 1;

d = mu./h(1:end-1) + mu./h(2:end);

d1 = -mu./h(2:end-1);

N = length(x)-2;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N, N);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

u = [0; A\b; 0];

u(find(x==0.54))


%% round 3
clc
clear
close all


mu = 1;

f = @(x) 2.*sin(x) + (1-x).*cos(x);

a = 0;
b = 2;

x = [0 0.33 0.57 0.63 0.88 1.07 1.31 1.4 1.79 1.91 2]';
h = diff(x);
% neumann a dx
BC = [4 0];

d = mu./h(1:end-1) + mu./h(2:end);

d = [d; mu/h(end)];

d1 = -mu./h(2:end);

N = length(x)-1;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N,N);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

b(1) = b(1) + mu./h(1) * BC(1);

b(end+1) = f(x(end))/2 * h(end)+BC(2);

u = [BC(1); A\b]
plot(x,u);

u(find(x==1.4))

%% 2

clc
clear
close all


mu = 1;
a = 0;
b = 1;

f = @(x) 3.*cos(x) + sqrt(x+2).*cos(x);

x = [0 0.13 0.2 0.32 0.42 0.58 0.67 0.75 0.8 0.94 1]';

h = diff(x);

d = mu./h(1:end-1) + mu./h(2:end);

d1 = -mu./h(2:end-1);

N = length(x)-2;
A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], N, N);

b = f(x(2:end-1))/2 .* (h(2:end) + h(1:end-1));

u = [0; A\b; 0]


plot(x,u)
u(find(x==0.42))


%% 3

clc
clear
close all

N = 30;

h = 1/N;

a = 0;
b = 1;
mu = 1;
gamma_f = @(x) (-2.*x).^2 * log(5)^2-2*log(5);
f = @(x) zeros(length(x),1);

% neumann dx e sx
BC = [0 -2/5*log(5)];

x = (a:h:b)';


uni = ones(length(x),1);
d = 2*mu./h .* uni;

d(1) = d(1)/2;
d(end) = d(end)/2;

d1 = -mu./h * uni;

A = spdiags([d, d1, d1], [0 1 -1], length(x), length(x));

xMedi = (x(2:end) + x(1:end-1))/2;

gammaVal = gamma_f(xMedi);

dR = h/3 * (gammaVal(1:end-1) + gammaVal(2:end));

dR = [h/3 * gammaVal(1); dR; h/3*gammaVal(end)];

dR1 = h/6 * gammaVal(1:end);

AR = spdiags([dR, [0;dR1], [dR1;0]], [0 1 -1], length(x), length(x));

A_tot = A+AR;

b = f(x(2:end-1))/2 .* (2*h);

b = [f(x(1))/2*h - BC(1); b; f(x(end))/2*h + BC(2)];

u = A_tot\b;

figure
hold on
plot(x,u)


u_esatta = @(x) 5.^(-x.^2);

plot(x, u_esatta(x))

zero = abs(u_esatta(0) - u(1))
uno = abs(u_esatta(1) - u(end))

%% 4
clc
clear
close all

mu = 1;
a = 0;
b = 1;

f = @(x) 2.*sin(x) + (1-x).*sin(x);

x = [0 0.1 0.29 0.35 0.47 0.58 0.69 0.77 0.88 0.93 1]';
h = diff(x);

d = mu./h(1:end-1) + mu./h(2:end);

d1 = -mu./h(2:end-1);

A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], length(x)-2, length(x)-2);

b = f(x(2:end-1))/2 .* (h(1:end-1) + h(2:end));

u = [0; A\b; 0];
plot(x,u)

u(find(x==0.69))









