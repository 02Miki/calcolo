%% 1

clc
clear
close all

f= @(x) -ones(length(x),1);
a = 0;
b = 1;

intervalli = 80;
h = (b-a)/intervalli;

x = a+h:h:b-h;

mu = 0.001;
N = length(x);
uni = ones(N,1);
d = (mu/h^2*2+1/h)*uni;

% ha un elemento in più che lascio per spdiags
d1 = -mu/h^2*uni;
d_1 = (-mu/h^2-1/h)*uni;

A = spdiags([d d1 d_1], [0 1 -1], N, N);

b = f(x);

u = A\b;

min(u)

% -0.9695, ok (CHECK APPUNTI, FILE SCARABOCCHI PER INFO)

%% 2

clc
clear
close all

f = @(x) cos(sin(x));

a = 0;
b = pi;

mu = 6;

BC = [1/2, 1/3];
intervalli = 1000;

h = (b-a)/intervalli;

x = a:h:b;

N = length(x)-2;
% -2 perché non voglio gli estremi, che già conosco per DD
uni = ones(N,1);

d = mu/h^2*2*uni;
d1 = mu/h^2*-uni;
d_1 = d1;


A = spdiags([d, d1, d_1], [0, 1, -1], N, N);

b = f(x(2:end-1));

b(1) = b(1) + mu/h^2 * BC(1);
b(end) = b(end) + mu/h^2 * BC(2);


u = [BC(1); A\b'; BC(2)];

max(u)

% 0.5699, ok

%% 3
clc
clear
close all


mu = 5;

% dirichlet omogeneo, le "ignoriamo"
BC = [0,0];

f = @(x) log(3.*x.^2+2);

intervalli = 100;
a = 0;
b = 1;

h = (b-a)/intervalli;

x = a+h:h:b-h;

N = length(x);

uni = ones(N, 1);
d = mu/h^2 * 2 * uni;
d1 = mu/h^2 * -1 * uni;

A = spdiags([d, d1, d1], [0 1 -1], N, N);

% non modifico gli estremi perché ho DD omogeneo
b = f(x);

% dovrei aggiungere anche le BC
u = A\b';

max(u)

% 0.0260, ok


%% 4
clc
clear
close all


a = 0;
b = 2*pi;
BC = [-1, 0];

mu = -1;

f = @(x) cos(x);

intervalli = 10;

h = (b-a)/intervalli;

% a+h perché il primo nodo lo so già
x = a+h:h:b;
N = length(x);
uni = ones(N-1, 1)

d = [mu/h^2 * 2 * uni; mu/h^2];
d1 = [mu/h^2*-uni; -mu/h^2];

A = spdiags([d, d1, d1], [0 1 -1], N, N);

b = f(x);
b(1) = b(1) + mu/h^2 * BC(1);


b(end) = 1/2*(b(end) + 2/h*BC(2));


u = [BC(1); A\b'];

plot([a, x], u)
u(find(x == pi) +1)


