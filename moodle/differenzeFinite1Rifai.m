%% 1

clc
clear
close all



a = 0;
b = 1;


BC = [0 0];


mu = 0.001;
f = @(x) -ones(length(x),1);

intervalli = 80;

h = (b-a)/intervalli;

d = mu/h^2 * (2) + 1/h;

d1 = mu/h^2 * (-1);

d_1 = mu/h^2 * (-1) + 1/h * (-1);

uni = ones(intervalli-1, 1)

A = spdiags(uni*[d, d1, d_1], [0 1 -1], intervalli-1, intervalli-1)

x = (a:h:b)';
b = f(x(2:end-1));

% non metto dirichlet perché è omogeneo

u = [BC(1); A\b; BC(2)];
min(u)



%% 2
clc
clear
close all

intervalli = 1000;

mu = 6;

f = @(x) cos(sin(x));

BC = [1/2 1/3];

a = 0;
b = pi;

h = (b-a)/intervalli;

d = mu/h^2 * 2;

d1 = mu/h^2 * (-1);


uni = ones(intervalli-1, 1);

A = spdiags(uni * [d, d1, d1], [0 1 -1], intervalli-1, intervalli-1);


x = (a:h:b)';
b = f(x(2:end-1));

b(1) = b(1) + mu/h^2 * BC(1);

b(end) = b(end) + mu/h^2 * BC(2);

u = [BC(1); A\b; BC(2)];

max(u)


%% 3

clc
clear
close all

f = @(x) log(3.*x.^2 + 2);

mu = 5;

BC = [0 0];

a = 0;
b = 1;
intervalli = 100;
h = (b-a)/intervalli;

d = mu/h^2 * (2);

d1 = mu/h^2 * (-1);

uni = ones(intervalli-1,1);

A = spdiags(uni *[d, d1, d1], [0 1 -1], intervalli-1, intervalli-1);

x = (a:h:b)';

b = f(x(2:end-1));


u = [BC(1); A\b; BC(2)];

max(u)

%% 4

clc
clear
close all

f = @(x) cos(x);

a = 0;
b = 2*pi;
intervalli = 10;
mu = -1;
h = (b-a)/intervalli;

% neumann a destra
BC = [-1 0];

% intervalli - 1 perché scrivo la roba di neumann extra "a mano"
uni = ones(intervalli-1, 1);

% aggiungo un altro elemento alla diagonale a dx per neumann, ma diviso 2
d = [mu/h^2 * (2) * uni; mu/h^2 * (1)];

% la diagonale è simmetrica, quindi d1 = d_1;
% aggiungo un 1 in più sempre per neumann
d1 = mu/h^2 * (-1) * [uni; 1];

A = spdiags([d, d1, d1], [0 1 -1], intervalli, intervalli);

x = (a:h:b)';

% da 2 a end perché non conosco end, quindi devo calcolarlo
b = f(x(2:end));

% dirichlet non omogeneo
b(1) = b(1) + mu/h^2 * BC(1);

% essendo BC(2) = 0 non importa, ma il meno è necessario perché nel testo
% dell'esercizio ci da u' = 0, mentre dovrebbe essere -u' = 0 per il caso
% "standard"
b(end) = b(end)/2 + BC(2)/h * (-1);

u = [BC(1); A\b];

plot(x, u)

u(find(x==pi))


%% extra

clc
clear
close all

f = @(x) sin(x);

a = 0;
b = 2*pi;
intervalli = 10;
mu = -1;
h = (b-a)/intervalli;

% neumann a destra
BC = [0 -1];

% intervalli - 1 perché scrivo la roba di neumann extra "a mano"
uni = ones(intervalli-1, 1);

% aggiungo un altro elemento alla diagonale a dx per neumann, ma diviso 2
d = [mu/h^2 * (2) * uni; mu/h^2 * (1)];

% la diagonale è simmetrica, quindi d1 = d_1;
% aggiungo un 1 in più sempre per neumann
d1 = mu/h^2 * (-1) * [uni; 1];

A = spdiags([d, d1, d1], [0 1 -1], intervalli, intervalli);

x = (a:h:b)';

% da 2 a end perché non conosco end, quindi devo calcolarlo
b = f(x(2:end));

% dirichlet non omogeneo
b(1) = b(1) + mu/h^2 * BC(1);

% essendo BC(2) = 0 non importa, ma il meno è necessario perché nel testo
% dell'esercizio ci da u' = 0, mentre dovrebbe essere -u' = 0 per il caso
% "standard"
b(end) = b(end)/2 + BC(2)/h * (-1);

u = [BC(1); A\b];

plot(x, u)

u(find(x==2*pi))


%% round 2

%% 1

clc
clear
close all

x0 = 0;
xf = 2*pi;

f = @(x) sin(x);
mu = -1;

% bisogna moltiplicare per mu perché psi è mu*u'. Essendo u' = -1, allora
% mu*u' = mu*-1
BC = [0 -1*mu];

intervalli = 10;

h = (xf-x0)/intervalli;

uni = ones(intervalli-1,1);

d = mu/h^2*2.*uni;

d = [d; mu/h^2];

d1 = -mu/h^2.*uni;

A = spdiags([d, [0;d1], [d1;0]], [0 1 -1], intervalli, intervalli);

x = (x0:h:xf)';
b = f(x(2:end));

b(1) = b(1) + mu/h^2*BC(1);
b(end) = b(end) + BC(2)/h;

u = [0; A\b];


u(end)

%% 2
clc
clear
close all

mu = 0.1;


--- da fare 


%% 3

clc
clear
close all

mu = 3;

f = @(x) exp(sin(x));

BC = [1 2];

a = -2;
b = 2;

intervalli = 1000;
h = (b-a)/intervalli;


uni = ones(intervalli-1,1);

d = 2*uni;
d1 = -1*uni;

A = spdiags(mu/h^2*[d, d1, d1], [0 1 -1], intervalli-1, intervalli-1);

x = (a:h:b)';
b = f(x(2:end-1));

b(1) = b(1) + mu/h^2 * BC(1);
b(end) = b(end) + mu/h^2 * BC(2);

u = [BC(1); A\b; BC(2)];

max(u)


%% 4

clc
clear
close all

mu = 2;

f = @(x) abs(cos(pi.*x));

a = 0;
b = 2;

BC = [0 0];
intervalli = 1000;

h = (b-a)/intervalli;

uni = ones(intervalli-1,1);
d = 2*uni;
d1 = -1*uni;

A = spdiags(mu/h^2.*[d, d1, d1], [0 1 -1], intervalli-1, intervalli-1);

x = (a:h:b)';
b = f(x(2:end-1));

u = [0; A\b; 0];

max(u)















