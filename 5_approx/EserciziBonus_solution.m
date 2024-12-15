%% Esercizio 1:

clc
clear
close all

x = 0:8; 
y = [4.9, 4.3, 7.1, 3.4, 2.9, 2.1, 3.5, 7.3, 2.3];

x_eval = 2.5;
approx = spline(x, [-0.5 y -2.0], x_eval) % Vanno aggiunte le condizioni 
% nel secondo input del comando spline

%% Esercizio 2:

clc
clear
close all

n = 250000;
e = ones(n,1);
A = spdiags([-1*e -1*e 4*e -1*e -1*e], [-500 -1 0 1 500], n,n);
b = ones(n, 1);
x = A \ b;

x(561)



%% Esercizio 3:

clc
clear
close all

f = @(x) exp(-1./(x+1));

n = 5;
t = - cos(0.5 * pi * ((2 * (0:n) + 1) / (n + 1)));

a = 0;
b = 1;

x = (a + b) * 0.5 + (b - a) * 0.5 * t;


c = polyfit(x, f(x), n) % Importante: vedasi l'help per 
% l'ordine in cui vengono ridati i coefficienti

%% Esercizio 4

clc
clear
close all

N = 20;
A = ones(N, N);
i = 1:N;
j = i';
A = A + exp(- abs(i - j));
b = zeros(N, 1);
b(1:2:N) = 1;

x = A \ b;
norm(x, inf)

%% Esercizio 5

clc
clear
close all

f = @(x) 1./(x.^4 + 1);

n = 6;
x = linspace(-4, 6, n+1);

c = polyfit(x, f(x), n);

x_eval = linspace(-4, 6, 500);
approx = polyval(c, x_eval);


figure
hold on
plot(x_eval, approx, 'r', 'linewidth', 2)
plot(x_eval, f(x_eval), 'b', 'linewidth', 2)
plot(x, f(x),'ro')


%% Esercizio 6

clc
clear
close all

% Dovete dare il path di dove è la vostra funzione di GaussSeidel, 
% oppure copiarlo e incollarlo nell'attuale cartella di lavoro
% addpath('/home/emma/Documenti/EMMA/Corsi/2024_02IHZMT_Metodi_numerici_e_calcolo_scientifico_274295/Lezione3')

n = 2500;
e = ones(n,1);
A = spdiags([-1*e -1*e 2*e -1*e -1*e], [-500 -1 0 1 500],n,n);
b = ones(n, 1);


x0 = zeros(n, 1);
tol = 1.0e-04;
kmax = 10000;
[x, res_rel, k] = myGaussSeidel(A, b, x0, tol, kmax);
k
M_GS = tril(A);
N_GS = A - M_GS;
rho_GS = myRho(M_GS,N_GS)


%% Esercizio 7

clc
clear
close all

A = [3.75 -1.25 -2.25 0.75;
    -1.25 3.75 0.75 -2.25;
    -2.25 0.75 3.75 -1.25;
    0.75 -2.25 -1.25 3.75];

[V,D] = eig(A);
V
diag(D)

w_1 = [5, -5, -5, -5]';
w_3 = [3, -3, 3, -3]';

w_3 / norm(w_3) % è una colonna di V
(A * w_3)./w_3 % corrisponde all'autovalore 2

%% Esercizio 8:

clc
clear
close all

xi = linspace(0, 2* pi, 21)';
f = @(x) x + sin(x).*cos(x) + 0.4 * (1.0 - sin(x).^(2));

V = [ones(size(xi)) xi cos(xi) sin(xi)];
y = f(xi);

% devo risolvere V*c = y

[Q,R] = qr(V);
c = R \ (Q' * y);

c